#include "TDMSSegment.h"

#include "TDMSExceptionHandler.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSSegment::TDMSSegment(std::ifstream& filestream, uint64_t currentSegment)
: m_FileStream(filestream)
, m_SegmentIndex(currentSegment)
, m_LeadIn(nullptr)
, m_RawDataPosition(filestream.tellg())
, m_NextSegmentPosition(filestream.tellg())
, m_NumberOfChunks(0)
, m_TotalSegmentDataSize(0)
{
  initializeSegment();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSSegment::~TDMSSegment() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSSegment::Pointer TDMSSegment::New(std::ifstream& filestream, uint64_t currentSegment)
{
  Pointer shared(new TDMSSegment(filestream, currentSegment));
  return shared;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSSegment::initializeSegment()
{
  m_LeadIn = TDMSLeadIn::New(m_FileStream);

  m_RawDataPosition += TDMSLeadIn::TDMSLEADINLENGTH + m_LeadIn->m_LeadInStruct.RawDataOffset;

  if(m_LeadIn->m_LeadInStruct.RemainingSegmentLength == 0xFFFFFFFFFFFFFFFF)
  {
    throw NonFatalTDMSException(TDMSExceptionMessages::EndOfFile);
  }

  m_NextSegmentPosition += TDMSLeadIn::TDMSLEADINLENGTH + m_LeadIn->m_LeadInStruct.RemainingSegmentLength;

  if(m_NextSegmentPosition < m_RawDataPosition)
  {
    std::string info("Next segment position (bytes): " + std::to_string(m_NextSegmentPosition) + "\n" + "Raw data position (bytes): " + std::to_string(m_RawDataPosition) + "\n" +
                     "Computed segment raw data size: " + std::to_string(m_NextSegmentPosition - m_RawDataPosition));
    throw FatalTDMSException(TDMSExceptionMessages::NegativeSegmentDataSize, info);
  }

  m_TotalSegmentDataSize = m_LeadIn->m_LeadInStruct.RemainingSegmentLength - m_LeadIn->m_LeadInStruct.RawDataOffset;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSSegment::readMetaData(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order)
{
  if(!m_LeadIn->m_ToCFlags.HasMetaData)
  {
    for(auto&& path : order)
    {
      objects[path]->m_MetaData->extendSegmentMetaData();
    }
    computeIncrementalChunks(objects, order);
    return;
  }

  if(!m_LeadIn->m_ToCFlags.HasNewObjList)
  {
    for(auto&& path : order)
    {
      objects[path]->m_MetaData->extendSegmentMetaData();
    }
  }

  char buffer[4];
  m_FileStream.read(buffer, 4);
  uint32_t numObjects = *(reinterpret_cast<uint32_t*>(buffer));

  std::list<std::string> objectsExtended;

  for(uint32_t i = 0; i < numObjects; i++)
  {
    m_FileStream.read(buffer, 4);
    uint32_t pathLength = *(reinterpret_cast<uint32_t*>(buffer));
    std::vector<char> pathBuffer(pathLength);
    m_FileStream.read(pathBuffer.data(), pathLength);
    std::string path(pathBuffer.data(), pathLength);

    TDMSObject::Pointer object = nullptr;
    if(objects.find(path) != objects.end())
    {
      object = objects[path];
      if(m_LeadIn->m_ToCFlags.HasNewObjList)
      {
        object->m_MetaData->extendSegmentMetaData();
      }
    }
    else
    {
      object = TDMSObject::New(path);
      object->m_MetaData->resizeSegmentMetaData(m_SegmentIndex + 1);
      objects[path] = object;
      order.push_back(path);
    }

    objectsExtended.push_back(path);
    object->m_MetaData->readSegmentMetaData(m_FileStream, m_SegmentIndex);
    object->populateMetaData();
  }

  if(m_LeadIn->m_ToCFlags.HasNewObjList)
  {
    while(objectsExtended.size() != objects.size())
    {
      auto iter = std::find_if_not(std::begin(objects), std::end(objects), [&objectsExtended](decltype(*std::begin(objects))& val) -> bool {
        auto pathIter = std::find(std::begin(objectsExtended), std::end(objectsExtended), val.first);
        return pathIter != objectsExtended.end();
      });
      if(iter == std::end(objects))
      {
        break;
      }
      iter->second->m_MetaData->extendSegmentMetaData(false);
      objectsExtended.push_back(iter->first);
    }
  }

  computeIncrementalChunks(objects, order);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSSegment::readRawData(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order)
{
  if(!m_LeadIn->m_ToCFlags.HasRawData)
  {
    return;
  }

  m_FileStream.seekg(m_RawDataPosition);

  for(uint64_t i = 0; i < m_NumberOfChunks; i++)
  {
    for(auto&& path : order)
    {
      if(objects[path]->m_MetaData->m_SegmentMetaData[m_SegmentIndex].HasData)
      {
        objects[path]->readRawData(m_FileStream, m_SegmentIndex);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSSegment::computeIncrementalChunks(std::unordered_map<std::string, TDMSObject::Pointer>& objects, std::vector<std::string>& order)
{
  uint64_t size = 0;

  for(auto&& path : order)
  {
    if(objects[path]->m_MetaData->m_SegmentMetaData[m_SegmentIndex].HasData)
    {
      size += objects[path]->m_MetaData->m_SegmentMetaData[m_SegmentIndex].TotalSegmentSize;
    }
  }

  if(size == 0)
  {
    if(m_TotalSegmentDataSize != size)
    {
      std::string info("Size indicated by meta data (bytes): " + std::to_string(size) + "\n" + "Size indicated by lead in (bytes): " + std::to_string(m_TotalSegmentDataSize));
      throw FatalTDMSException(TDMSExceptionMessages::SegmentDataSizeMismatch, info);
    }
    m_NumberOfChunks = 0;
    return;
  }

  if((m_TotalSegmentDataSize % size) != 0)
  {
    std::string info("Size indicated by meta data (bytes): " + std::to_string(size) + "\n" + "Size indicated by lead in (bytes): " + std::to_string(m_TotalSegmentDataSize) + "\n" +
                     "Computed number of chunks: " + std::to_string(static_cast<float>(m_TotalSegmentDataSize) / static_cast<float>(size)));
    throw FatalTDMSException(TDMSExceptionMessages::SegmentChunkMismatch, info);
  }

  m_NumberOfChunks = m_TotalSegmentDataSize / size;

  for(auto&& path : order)
  {
    if(objects[path]->m_MetaData->m_SegmentMetaData[m_SegmentIndex].HasData)
    {
      objects[path]->updateSizeInformation(m_NumberOfChunks, m_SegmentIndex);
    }
  }
}
