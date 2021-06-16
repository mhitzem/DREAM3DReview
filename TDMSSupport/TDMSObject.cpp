#include "TDMSObject.h"

#include "TDMSExceptionHandler.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSObject::TDMSObject(const std::string& path)
: m_Path(path)
, m_MetaData(TDMSMetaData::New())
, m_DataType(nullptr)
, m_ArrayDimension(0)
, m_NumberOfValues(0)
, m_TotalSize(0)
, m_CurrentDataPosition(0)
, m_HasData(false)
, m_HasInitializedMetaData(false)
, m_Data(nullptr)
{
  determineObjectType();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSObject::~TDMSObject() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSObject::Pointer TDMSObject::New(const std::string& path)
{
  Pointer shared(new TDMSObject(path));
  return shared;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::populateMetaData()
{
  if(m_MetaData->m_SegmentMetaData.empty())
  {
    return;
  }
  if(!m_HasInitializedMetaData)
  {
    m_DataType = m_MetaData->m_SegmentMetaData.back().DataType;
    m_ArrayDimension = m_MetaData->m_SegmentMetaData.back().ArrayDimension;
    m_HasData = m_MetaData->m_SegmentMetaData.back().HasData;
    m_HasInitializedMetaData = true;
  }
  else
  {
    if(m_DataType->name() != m_MetaData->m_SegmentMetaData.back().DataType->name())
    {
      std::string info("Data type stored in object: " + m_DataType->name() + "\n" + "Data type indicated by meta data: " + m_MetaData->m_SegmentMetaData.back().DataType->name());
      throw FatalTDMSException(TDMSExceptionMessages::ObjectMetaDataMismatch, info);
    }
    if(m_ArrayDimension != m_MetaData->m_SegmentMetaData.back().ArrayDimension)
    {
      std::string info("Array dimension stored in object: " + std::to_string(m_ArrayDimension) + "\n" +
                       "Array dimension indicated by meta data: " + std::to_string(m_MetaData->m_SegmentMetaData.back().ArrayDimension));
      throw FatalTDMSException(TDMSExceptionMessages::ObjectMetaDataMismatch, info);
    }
    if(m_MetaData->m_SegmentMetaData.back().HasData)
    {
      m_HasData = true;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::updateSizeInformation(uint64_t numChunks, size_t index)
{
  m_NumberOfValues += (m_MetaData->m_SegmentMetaData[index].NumberOfValues * numChunks);
  m_TotalSize = m_NumberOfValues * m_DataType->size() * m_ArrayDimension;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::generateDataArray()
{
  if(m_DataType && m_HasData)
  {
    m_Data = m_DataType->generateDataArray(m_NumberOfValues, m_BaseName);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::allocate()
{
  if(m_Data && m_DataType && m_HasData)
  {
    m_DataType->allocateDataArray(m_Data);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::readRawData(std::ifstream& filestream, uint64_t index)
{
  if(m_Data && m_DataType && m_HasData)
  {
    m_DataType->readArrayFromFile(filestream, m_Data, m_CurrentDataPosition, m_MetaData->m_SegmentMetaData[index].NumberOfValues);
    m_CurrentDataPosition += m_MetaData->m_SegmentMetaData[index].NumberOfValues;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::string TDMSObject::parseChannelName()
{
  std::string::size_type pos = m_Path.rfind("'", m_Path.size() - 2);
  std::string channel = m_Path.substr(pos + 1, m_Path.size() - pos - 2);
  return channel;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::string TDMSObject::parseGroupName()
{
  std::string group = m_Path.substr(2);
  group.pop_back();
  return group;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSObject::determineObjectType()
{
  if(m_Path == "/")
  {
    m_ObjectType = Type::Root;
    m_BaseName = m_Path;
  }
  else if(std::count(std::begin(m_Path), std::end(m_Path), '/') == 1)
  {
    m_ObjectType = Type::Group;
    m_BaseName = parseGroupName();
  }
  else
  {
    m_ObjectType = Type::Channel;
    m_BaseName = parseChannelName();
  }
}