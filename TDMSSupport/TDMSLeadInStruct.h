#ifndef _tdmsleadinstruct_h_
#define _tdmsleadinstruct_h_

#include <cstdint>

#pragma pack(push)
#pragma pack(1)

struct TDMSLeadInStruct
{
  char TDMSTag[4];
  uint32_t TableOfContents;
  uint32_t VersionNumber;
  uint64_t RemainingSegmentLength;
  uint64_t RawDataOffset;
};

#pragma pack(pop)

#endif
