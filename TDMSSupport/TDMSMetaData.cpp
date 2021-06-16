#include "TDMSMetaData.h"

#include "TDMSDataTypeFactory.h"
#include "TDMSExceptionHandler.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSMetaData::TDMSMetaData()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSMetaData::~TDMSMetaData() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSMetaData::Pointer TDMSMetaData::New()
{
  Pointer shared(new TDMSMetaData());
  return shared;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSMetaData::extendSegmentMetaData()
{
  if(m_SegmentMetaData.empty())
  {
    m_SegmentMetaData.resize(1);
  }
  else
  {
    m_SegmentMetaData.resize(m_SegmentMetaData.size() + 1, m_SegmentMetaData.back());
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSMetaData::extendSegmentMetaData(size_t count)
{
  m_SegmentMetaData.resize(m_SegmentMetaData.size() + count);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSMetaData::extendSegmentMetaData(bool hasData)
{
  extendSegmentMetaData();
  m_SegmentMetaData.back().HasData = hasData;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSMetaData::resizeSegmentMetaData(size_t size)
{
  m_SegmentMetaData.resize(size);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSMetaData::readSegmentMetaData(std::ifstream& filestream, size_t index)
{
  char buffer[8];
  filestream.read(buffer, 4);
  m_SegmentMetaData[index].RawDataIndex = *(reinterpret_cast<uint32_t*>(buffer));

  if(m_SegmentMetaData[index].RawDataIndex == 0xFFFFFFFF)
  {
    m_SegmentMetaData[index].HasData = false;
  }
  else if(m_SegmentMetaData[index].RawDataIndex == 0x00000000)
  {
    m_SegmentMetaData[index].HasData = true;
  }
  else
  {
    m_SegmentMetaData[index].HasData = true;

    filestream.read(buffer, 4);
    uint32_t dataTypeIndex = *(reinterpret_cast<uint32_t*>(buffer));

    TDMSDataTypeFactory* factory = TDMSDataTypeFactory::Instance();
    TDMSDataType::Pointer dataType = factory->getDataType(dataTypeIndex);
    m_SegmentMetaData[index].DataType = dataType;

    filestream.read(buffer, 4);
    m_SegmentMetaData[index].ArrayDimension = *(reinterpret_cast<uint32_t*>(buffer));
    if(m_SegmentMetaData[index].ArrayDimension != 1)
    {
      std::string info("Array dimension from TDMS file: " + std::to_string(m_SegmentMetaData[index].ArrayDimension) + "\n" + "Required array dimension: 1");
      throw FatalTDMSException(TDMSExceptionMessages::UnexpectedArrayDimension, info);
    }

    filestream.read(buffer, 8);
    m_SegmentMetaData[index].NumberOfValues = *(reinterpret_cast<uint64_t*>(buffer));

    if(m_SegmentMetaData[index].DataType->name() == "tdsTypeString")
    {
      filestream.read(buffer, 8);
      m_SegmentMetaData[index].TotalSegmentSize = *(reinterpret_cast<uint64_t*>(buffer));
    }
    else
    {
      m_SegmentMetaData[index].TotalSegmentSize = m_SegmentMetaData[index].NumberOfValues * dataType->size() * m_SegmentMetaData[index].ArrayDimension;
    }
  }

  filestream.read(buffer, 4);
  uint32_t numProperties = *(reinterpret_cast<uint32_t*>(buffer));

  for(uint32_t i = 0; i < numProperties; i++)
  {
    filestream.read(buffer, 4);
    uint32_t nameLength = *(reinterpret_cast<uint32_t*>(buffer));
    std::vector<char> nameBuffer(nameLength);
    filestream.read(nameBuffer.data(), nameLength);
    std::string name(nameBuffer.data(), nameLength);

    filestream.read(buffer, 4);
    uint32_t dataTypeIndex = *(reinterpret_cast<uint32_t*>(buffer));

    TDMSDataTypeFactory* factory = TDMSDataTypeFactory::Instance();
    TDMSDataType::Pointer dataType = factory->getDataType(dataTypeIndex);

    TDMSProperty::Pointer property = TDMSProperty::New(dataType, dataType->readSingleValueFromFile(filestream, name));
    m_Properties[name] = property;
  }
}