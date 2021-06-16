#include "TDMSFileProxy.h"

#include "TDMSExceptionHandler.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSFileProxy::TDMSFileProxy(const std::string& file)
: m_File(file)
, m_FileStream(std::ifstream(m_File.data(), std::ios::binary | std::ios::in))
, m_ObjectsAllocated(false)
, m_MetaDataRead(false)
{
  if(!m_FileStream.good())
  {
    throw FatalTDMSException(TDMSExceptionMessages::BadFile);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSFileProxy::~TDMSFileProxy() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSFileProxy::Pointer TDMSFileProxy::New(const std::string& file)
{
  Pointer shared(new TDMSFileProxy(file));
  return shared;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSFileProxy::readMetaData()
{
  if(m_MetaDataRead)
  {
    return;
  }
  uint64_t currentSegment = 0;
  while(true)
  {
    TDMSSegment::Pointer segment = nullptr;

    try
    {
      segment = TDMSSegment::New(m_FileStream, currentSegment);
    } catch(const NonFatalTDMSException&)
    {
      break;
    }
    if(!segment)
    {
      break;
    }

    segment->readMetaData(m_Objects, m_ObjectOrder);
    m_Segments.push_back(segment);
    currentSegment++;
    if(m_FileStream.eof())
    {
      break;
    }
    m_FileStream.seekg(segment->m_NextSegmentPosition);
  }

  m_FileStream.clear();

  for(auto&& path : m_ObjectOrder)
  {
    m_Objects[path]->generateDataArray();
  }
  m_MetaDataRead = true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSFileProxy::readRawData()
{
  if(!m_MetaDataRead)
  {
    readMetaData();
    allocateObjects();
  }
  for(auto&& segment : m_Segments)
  {
    segment->readRawData(m_Objects, m_ObjectOrder);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSFileProxy::allocateObjects()
{
  if(m_ObjectsAllocated)
  {
    return;
  }
  if(!m_MetaDataRead)
  {
    readMetaData();
  }
  for(auto&& path : m_ObjectOrder)
  {
    m_Objects[path]->allocate();
  }
  m_ObjectsAllocated = true;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSObject::Pointer TDMSFileProxy::rootObject()
{
  return (std::begin(extractObjectsOfType(TDMSObject::Type::Root))->second);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::unordered_map<std::string, TDMSObject::Pointer> TDMSFileProxy::groupObjects()
{
  return extractObjectsOfType(TDMSObject::Type::Group);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::unordered_map<std::string, TDMSObject::Pointer> TDMSFileProxy::channelObjects()
{
  return extractObjectsOfType(TDMSObject::Type::Channel);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::unordered_map<std::string, TDMSObject::Pointer> TDMSFileProxy::extractObjectsOfType(TDMSObject::Type type)
{
  std::unordered_map<std::string, TDMSObject::Pointer> objects;

  std::for_each(std::begin(m_Objects), std::end(m_Objects), [&objects, &type](decltype(*std::begin(m_Objects))& val) {
    if(val.second->objectType() == type)
    {
      objects[val.second->baseName()] = val.second;
    }
  });

  return objects;
}
