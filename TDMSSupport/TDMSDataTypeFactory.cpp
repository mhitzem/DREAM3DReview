#include "TDMSDataTypeFactory.h"

#include "TDMSExceptionHandler.h"

TDMSDataTypeFactory* TDMSDataTypeFactory::self = nullptr;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSDataTypeFactory::TDMSDataTypeFactory()
{
  initializeDataTypes();
  self = this;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSDataTypeFactory::~TDMSDataTypeFactory() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSDataTypeFactory* TDMSDataTypeFactory::Instance()
{
  if(self == nullptr)
  {
    self = new TDMSDataTypeFactory();
  }
  return self;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSDataType::Pointer TDMSDataTypeFactory::getDataType(uint32_t index)
{
  TDMSDataType::Pointer dataType = nullptr;
  try
  {
    dataType = m_DataTypes.at(index);
  } catch(const std::out_of_range&)
  {
    std::string info("Data type index from file: " + std::to_string(index) + "\n" + "Supported data type indices: \n" + getListOfSupportedDataTypes());
    throw FatalTDMSException(TDMSExceptionMessages::UnsupportedDataType, info);
  }
  return dataType;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSDataTypeFactory::initializeDataTypes()
{
  m_DataTypes[1] = TDMSDataType::New("tdsTypeI8", 1, TDMSDataTypeHelpers::ValueReaderFactory<int8_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<int8_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<int8_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<int8_t>());
  m_DataTypes[2] = TDMSDataType::New("tdsTypeI16", 2, TDMSDataTypeHelpers::ValueReaderFactory<int16_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<int16_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<int16_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<int16_t>());
  m_DataTypes[3] = TDMSDataType::New("tdsTypeI32", 4, TDMSDataTypeHelpers::ValueReaderFactory<int32_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<int32_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<int32_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<int32_t>());
  m_DataTypes[4] = TDMSDataType::New("tdsTypeI64", 8, TDMSDataTypeHelpers::ValueReaderFactory<int64_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<int64_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<int64_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<int64_t>());
  m_DataTypes[5] = TDMSDataType::New("tdsTypeU8", 1, TDMSDataTypeHelpers::ValueReaderFactory<uint8_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<uint8_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<uint8_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<uint8_t>());
  m_DataTypes[6] = TDMSDataType::New("tdsTypeU16", 2, TDMSDataTypeHelpers::ValueReaderFactory<uint16_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<uint16_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<uint16_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<uint16_t>());
  m_DataTypes[7] = TDMSDataType::New("tdsTypeU32", 4, TDMSDataTypeHelpers::ValueReaderFactory<uint32_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<uint32_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<uint32_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<uint32_t>());
  m_DataTypes[8] = TDMSDataType::New("tdsTypeU64", 8, TDMSDataTypeHelpers::ValueReaderFactory<uint64_t>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<uint64_t>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<uint64_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<uint64_t>());
  m_DataTypes[9] = TDMSDataType::New("tdsTypeSingleFloat", 4, TDMSDataTypeHelpers::ValueReaderFactory<float>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<float>(),
                                     TDMSDataTypeHelpers::ArrayAllocatorFactory<float>(), TDMSDataTypeHelpers::ArrayReaderFactory<float>());
  m_DataTypes[10] = TDMSDataType::New("tdsTypeDoubleFloat", 8, TDMSDataTypeHelpers::ValueReaderFactory<double>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<double>(),
                                      TDMSDataTypeHelpers::ArrayAllocatorFactory<double>(), TDMSDataTypeHelpers::ArrayReaderFactory<double>());
  m_DataTypes[0x20] = TDMSDataType::New("tdsTypeString", 0, TDMSDataTypeHelpers::StringReaderFactory(), TDMSDataTypeHelpers::StringArrayGeneratorFactory(),
                                        TDMSDataTypeHelpers::StringArrayAllocatorFactory(), TDMSDataTypeHelpers::StringArrayReaderFactory());
  m_DataTypes[0x21] = TDMSDataType::New("tdsTypeBoolean", 1, TDMSDataTypeHelpers::ValueReaderFactory<bool>(), TDMSDataTypeHelpers::ArrayGeneratorFactory<bool>(),
                                        TDMSDataTypeHelpers::ArrayAllocatorFactory<bool>(), TDMSDataTypeHelpers::ArrayReaderFactory<bool>());
  m_DataTypes[0x44] = TDMSDataType::New("tdsTypeTimeStamp", 16, TDMSDataTypeHelpers::TimeStampReaderFactory(), TDMSDataTypeHelpers::ArrayGeneratorFactory<uint8_t>(),
                                        TDMSDataTypeHelpers::ArrayAllocatorFactory<uint8_t>(), TDMSDataTypeHelpers::ArrayReaderFactory<uint8_t>());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::string TDMSDataTypeFactory::getListOfSupportedDataTypes()
{
  std::string types;

  for(auto&& pair : m_DataTypes)
  {
    types += (std::to_string(pair.first) + " (" + pair.second->name() + ")\n");
  }

  return types;
}