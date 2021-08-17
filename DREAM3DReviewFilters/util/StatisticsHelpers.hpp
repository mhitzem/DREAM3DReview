#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

namespace StatisticsHelpers
{
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
T findMin(const C<T, Ts...>& source)
{
  if(source.empty())
  {
    return static_cast<T>(0);
  }
  return (*std::min_element(std::cbegin(source), std::cend(source)));
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
T findMax(const C<T, Ts...>& source)
{
  if(source.empty())
  {
    return static_cast<T>(0);
  }
  return (*std::max_element(std::cbegin(source), std::cend(source)));
}

// -----------------------------------------------------------------------------
template <class Container>
auto computeSum(const Container& source)
{
  using T = typename Container::value_type;
  if constexpr(std::is_integral_v<T>)
  {
    if constexpr(std::is_signed_v<T>)
    {
      return std::accumulate(std::cbegin(source), std::cend(source), static_cast<int64_t>(0));
    }
    else
    {
      return std::accumulate(std::cbegin(source), std::cend(source), static_cast<uint64_t>(0));
    }
  }
  else // Use Kahan Summation for float/double
  {
    double sum = 0.0;

    // Variable to store the error
    double c = 0.0;

    // Loop to iterate over the array
    for(const T& f : source)
    {
      double y = static_cast<double>(f) - c;
      double t = sum + y;

      // Algebraically, c is always 0
      // when t is replaced by its
      // value from the above expression.
      // But, when there is a loss,
      // the higher-order y is cancelled
      // out by subtracting y from c and
      // all that remains is the
      // lower-order error in c
      c = (t - sum) - y;
      sum = t;
    }
    return sum;
  }
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findMean(const C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = static_cast<float>(computeSum(source));

  return sum / static_cast<float>(source.size());
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts>
bool findMean(const C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::cbegin(source), std::cend(source), true);
  return count >= (source.size() - count);
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findMedian(const C<T, Ts...>& source)
{
  // Need a copy, not a reference, since we will be sorting the input array
  std::vector<T> tmpList{std::cbegin(source), std::cend(source)};
  if(tmpList.empty())
  {
    return 0.0f;
  }
  std::sort(tmpList.begin(), tmpList.end());
  float medVal = 0.0f;
  if(tmpList.size() % 2 == 1)
  {
    size_t halfElements = static_cast<size_t>(std::floor(tmpList.size() / 2.0f));
    medVal = tmpList[halfElements];
  }
  else
  {
    size_t idxLow = (tmpList.size() / 2) - 1;
    size_t idxHigh = tmpList.size() / 2;
    medVal = (tmpList[idxLow] + tmpList[idxHigh]) * 0.5f;
  }
  return medVal;
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
float findStdDeviation(const C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  std::vector<double> difference(source.size());
  float sum = static_cast<float>(computeSum(source));
  float mean = static_cast<double>(sum / source.size());
  std::transform(std::cbegin(source), std::cend(source), std::begin(difference), [mean](float a) { return a - mean; });
  float squaredSum = std::inner_product(std::cbegin(difference), std::cend(difference), std::cbegin(difference), 0.0f);
  return std::sqrt(squaredSum / source.size());
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts>
bool findStdDeviation(const C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::cbegin(source), std::cend(source), true);
  return count >= (source.size() - count);
}

// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts>
double findSummation(const C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = static_cast<float>(computeSum(source));
  return sum;
}
} // namespace StatisticsHelpers

#ifdef STATISTICS_FILTER_CLASS_NAME
namespace
{
// Create Enumerations to allow the created Attribute Arrays to take part in renaming
enum createdPathID : RenameDataPath::DataID_t
{
  DataArrayID30 = 30, // Length
  DataArrayID31 = 31, // Min
  DataArrayID32 = 32, // Max
  DataArrayID33 = 33, // Mean
  DataArrayID34 = 34, // Median
  DataArrayID35 = 35, // StdDev
  DataArrayID36 = 36, // Summation
  DataArrayID37 = 37, // Histogram
  DataArrayID38 = 38, //
  DataArrayID39 = 39, // StandardizedArray
  DataArrayID40 = 40, //
};
} // namespace

/**
 * @brief createCompatibleArrays Creates the output statistics arrays with compatible types based on the
 * user options and incoming DataArray type
 */
template <typename T>
void STATISTICS_FILTER_CLASS_NAME::createCompatibleArrays(QVector<DataArrayPath>& dataArrayPaths)
{
  std::vector<size_t> cDims = {1};
  using DataArrayType = DataArray<T>;

  if(m_FindLength)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getLengthArrayName());
    m_LengthPtr = getDataContainerArray()->createNonPrereqArrayFromPath<MeshIndexArrayType>(this, path, 0, cDims, "", DataArrayID30);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMin)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMinimumArrayName());
    m_MinimumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArrayType>(this, path, 0, cDims, "", DataArrayID31);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMax)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMaximumArrayName());
    m_MaximumPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArrayType>(this, path, 0, cDims, "", DataArrayID32);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMean)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMeanArrayName());
    m_MeanPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID33);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindMedian)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getMedianArrayName());
    m_MedianPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID34);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindStdDeviation)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getStdDeviationArrayName());
    m_StandardDeviationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<FloatArrayType>(this, path, 0, cDims, "", DataArrayID35);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  if(m_FindSummation)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getSummationArrayName());
    m_SummationPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>>(this, path, 0, cDims, "", DataArrayID36);
    if(getErrorCode() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }
}
#endif
