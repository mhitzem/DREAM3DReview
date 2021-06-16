#include "TDMSLeadIn.h"

#include <string>

#include "TDMSExceptionHandler.h"

const std::string TDMSLeadIn::TDMSTAG = "TDSm";
const std::set<uint32_t> TDMSLeadIn::TDMSVERSIONS = {4712, 4713};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSLeadIn::TDMSLeadIn(std::ifstream& filestream)
: m_FileStream(filestream)
{
  constructLeadIn();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSLeadIn::~TDMSLeadIn() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TDMSLeadIn::Pointer TDMSLeadIn::New(std::ifstream& filestream)
{
  Pointer shared(new TDMSLeadIn(filestream));
  return shared;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TDMSLeadIn::constructLeadIn()
{
  char buffer[TDMSLEADINLENGTH];
  m_FileStream.read(buffer, TDMSLEADINLENGTH);
  if(m_FileStream.eof())
  {
    throw NonFatalTDMSException(TDMSExceptionMessages::EndOfFile);
  }
  m_LeadInStruct = *(reinterpret_cast<TDMSLeadInStruct*>(buffer));

  std::string tdmstag(m_LeadInStruct.TDMSTag, 4);
  if(tdmstag != TDMSTAG)
  {
    std::string info("Tag from TDMS file: " + tdmstag + "\n" + "Valid TDMS tag: " + TDMSTAG);
    throw FatalTDMSException(TDMSExceptionMessages::InvalidTag, info);
  }

  if(TDMSVERSIONS.find(m_LeadInStruct.VersionNumber) == TDMSVERSIONS.end())
  {
    std::string info("Version number from TDMS file: " + std::to_string(m_LeadInStruct.VersionNumber) + "\n" + "Valid TDMS version numbers: 4712; 4713");
    throw FatalTDMSException(TDMSExceptionMessages::InvalidVersion, info);
  }

  ToCBitMasks masks;
  m_ToCFlags.HasMetaData = (m_LeadInStruct.TableOfContents & masks.MetaData) != 0;
  m_ToCFlags.HasRawData = (m_LeadInStruct.TableOfContents & masks.RawData) != 0;
  m_ToCFlags.HasDAQmxRawData = (m_LeadInStruct.TableOfContents & masks.DAQmxRawData) != 0;
  m_ToCFlags.IsInterleavedData = (m_LeadInStruct.TableOfContents & masks.InterleavedData) != 0;
  m_ToCFlags.IsBigEndian = (m_LeadInStruct.TableOfContents & masks.BigEndian) != 0;
  m_ToCFlags.HasNewObjList = (m_LeadInStruct.TableOfContents & masks.NewObjList) != 0;

  if(m_ToCFlags.IsBigEndian)
  {
    throw FatalTDMSException(TDMSExceptionMessages::IsBigEndian);
  }
  if(m_ToCFlags.IsInterleavedData)
  {
    throw FatalTDMSException(TDMSExceptionMessages::IsInterleaved);
  }
  if(m_ToCFlags.HasDAQmxRawData)
  {
    throw FatalTDMSException(TDMSExceptionMessages::HasDAQmx);
  }
}
