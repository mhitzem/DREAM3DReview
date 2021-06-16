#ifndef _tdmsexceptionhandler_h_
#define _tdmsexceptionhandler_h_

namespace TDMSExceptionMessages
{
const std::string TDMSFileError("TDMS File Error: ");
const std::string TDMSLeadInError("TDMS Lead In Error: ");
const std::string TDMSMetaDataError("TDMS Meta Data Error: ");
const std::string TDMSDataTypeError("TDMS Data Type Error: ");
const std::string BadFile = TDMSFileError + "Unable to open file";
const std::string EndOfFile = TDMSFileError + "Reached end of file";
const std::string InvalidTag = TDMSLeadInError + "Lead in contains invalid tag";
const std::string InvalidVersion = TDMSLeadInError + "Lead in contains invalid version number";
const std::string IsBigEndian = TDMSLeadInError + "Lead in indicates data are big endian; only little endian data are supported";
const std::string IsInterleaved = TDMSLeadInError + "Lead in indicates data are interleaved; only contiguous data are supported";
const std::string HasDAQmx = TDMSLeadInError + "Lead in indicates segment contains DAQmx data; importing DAQmx data is not supported";
const std::string NegativeSegmentDataSize = TDMSLeadInError + "Lead in byte offset values indicate negative data size for segment";
const std::string SegmentDataSizeMismatch = TDMSMetaDataError + "Lead in indicates segment contains raw data, but no objects in the segment have meta data with associated raw data";
const std::string SegmentChunkMismatch = TDMSMetaDataError + "Meta data for objects in the segment indicates a number of objects that is not an integer factor of the total segment raw data size";
const std::string UnexpectedArrayDimension = TDMSMetaDataError + "Meta data for object indicates array dimension other than 1; only scalar dimension arrays are supported";
const std::string UnsupportedDataType = TDMSDataTypeError + "Encountered an unsupported TDMS data type";
const std::string ObjectMetaDataMismatch = TDMSMetaDataError + "Meta data for same object in multiple segments does not contain matching data type or array dimension";
} // namespace TDMSExceptionMessages

class NonFatalTDMSException : public std::exception
{
public:
  NonFatalTDMSException(const std::string& message = "Undefined non-fatal TDMS exception", const std::string& info = "")
  : std::exception()
  , m_Message(message)
  , m_Information(info)
  {
    constructMessageWithInformation();
  }

  virtual ~NonFatalTDMSException() = default;

  NonFatalTDMSException(const NonFatalTDMSException& exc)
  {
    m_Message = exc.getMessage();
    m_Information = exc.getInformation();
    constructMessageWithInformation();
  }

  std::string getMessage() const
  {
    return m_Message;
  }

  std::string getInformation() const
  {
    return m_Information;
  }

  const char* what() const noexcept
  {
    static std::string message = m_MessageWithInformation;
    return message.c_str();
  }

private:
  void constructMessageWithInformation()
  {
    m_MessageWithInformation = m_Message;
    if(m_Information.length() == 0)
    {
      return;
    }

    m_MessageWithInformation += ("\n" + m_Information);
  }

  std::string m_Message;
  std::string m_Information;
  std::string m_MessageWithInformation;
};

class FatalTDMSException : public std::exception
{
public:
  FatalTDMSException(const std::string& message = "Undefined fatal TDMS exception", const std::string& info = "")
  : std::exception()
  , m_Message(message)
  , m_Information(info)
  {
    constructMessageWithInformation();
  }

  virtual ~FatalTDMSException() = default;

  FatalTDMSException(const FatalTDMSException& exc)
  {
    m_Message = exc.getMessage();
    m_Information = exc.getInformation();
    constructMessageWithInformation();
  }

  std::string getMessage() const
  {
    return m_Message;
  }

  std::string getInformation() const
  {
    return m_Information;
  }

  const char* what() const noexcept
  {
    static std::string message = m_MessageWithInformation;
    return message.c_str();
  }

private:
  void constructMessageWithInformation()
  {
    m_MessageWithInformation = m_Message;
    if(m_Information.length() == 0)
    {
      return;
    }

    m_MessageWithInformation += ("\n" + m_Information);
  }

  std::string m_Message;
  std::string m_Information;
  std::string m_MessageWithInformation;
};

#endif /* _tdmsexceptionhandler_h_ */
