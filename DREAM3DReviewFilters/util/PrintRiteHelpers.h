#pragma once

#include <cassert>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Dense>

#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/TDMSSupport/TDMSObject.h"

namespace PrintRiteHelpers
{
class PrintRiteChannels
{
public:
  PrintRiteChannels()
  {
    m_HfChannelNames.insert("HF Time");
    m_HfChannelNames.insert("HF Temperature");
    m_HfChannelNames.insert("Laser Drive");
    m_HfChannelNames.insert("Audio");
    m_HfChannelNames.insert("Photodiode");
    m_HfChannelNames.insert("On-Axis PD");
    m_HfChannelNames.insert("Off-Axis PD");
    m_HfChannelNames.insert("X Position");
    m_HfChannelNames.insert("Y Position");
    m_HfChannelNames.insert("Laser Power");
    m_LfChannelNames.insert("LF Time");
    m_LfChannelNames.insert("LF Temperature");
    m_LfChannelNames.insert("O2");
  }

  enum class ChannelType
  {
    HF,
    LF,
    Unknown
  };

  std::unordered_set<std::string> getHfArrayNames()
  {
    return m_HfChannelNames;
  }

  std::unordered_set<std::string> getLfArrayNames()
  {
    return m_LfChannelNames;
  }

  std::vector<BoolArrayType::Pointer> thresholdHfChannels(std::unordered_map<std::string, TDMSObject::Pointer> channelObjects, const std::vector<std::string>& channels, const double& tolerance)
  {
    size_t numValidObjects = 0;
    std::vector<std::string> foundHfArrayNames;
    for(auto&& name : m_HfChannelNames)
    {
      if(std::find(std::begin(channels), std::end(channels), name) != std::end(channels) && channelObjects.count(name) > 0)
      {
        numValidObjects++;
        foundHfArrayNames.push_back(name);
      }
    }

    std::vector<typename BoolArrayType::Pointer> thresholdArrays(numValidObjects);
    for(auto i = 0; i < foundHfArrayNames.size(); i++)
    {
      TDMSObject::Pointer tdmsObject = channelObjects[foundHfArrayNames[i]];
      IDataArray::Pointer data = tdmsObject->data();
      // TODO: The threshold array name should be an option, but for now hard set it since we
      //      only plan on using thresholding the laser drive
      BoolArrayType::Pointer thresholdData = BoolArrayType::CreateArray(data->getNumberOfTuples(), std::string("Laser On"), true);
      thresholdData->initializeWithValue(false);
      bool* thresholdDataPtr = thresholdData->getPointer(0);
      TDMSDataType::Pointer dataType = tdmsObject->dataType();
      // TODO: For now, assert that all types are double; eventually this can be made generic for all types,
      //      but we know that the data coming from the PrintRite are all stored as doubles
      assert(dataType->name() == "tdsTypeDoubleFloat");
      DoubleArrayType::Pointer tdmsData = std::dynamic_pointer_cast<DoubleArrayType>(data);
      double* tdmsDataPtr = tdmsData->getPointer(0);
      for(size_t j = 0; j < data->getNumberOfTuples(); j++)
      {
        thresholdDataPtr[j] = (tdmsDataPtr[j] > tolerance);
      }
      thresholdArrays[i] = thresholdData;
    }
    return thresholdArrays;
  }

  std::vector<IDataArray::Pointer> getChannelsOfType(ChannelType type, std::unordered_map<std::string, TDMSObject::Pointer> channelObjects)
  {
    std::vector<std::string> validChannels = matchChannelsOfType(type, channelObjects, std::vector<std::string>());
    std::vector<IDataArray::Pointer> channels;
    for(auto&& pair : channelObjects)
    {
      if(std::find(std::begin(validChannels), std::end(validChannels), pair.first) != std::end(validChannels))
      {
        channels.push_back(pair.second->data());
      }
    }
    return channels;
  }

  template <typename T>
  std::vector<typename DataArray<T>::Pointer> castChannelsTo(ChannelType type, std::unordered_map<std::string, TDMSObject::Pointer> channelObjects,
                                                             const std::vector<std::string>& ignoreList = std::vector<std::string>())
  {
    // size_t numValidObjects = 0;
    std::vector<std::string> validChannels = matchChannelsOfType(type, channelObjects, ignoreList);
    std::vector<typename DataArray<T>::Pointer> downcastArrays(validChannels.size());
    for(auto i = 0; i < validChannels.size(); i++)
    {
      TDMSObject::Pointer tdmsObject = channelObjects[validChannels[i]];
      IDataArray::Pointer data = tdmsObject->data();
      typename DataArray<T>::Pointer downcastData = DataArray<T>::CreateArray(data->getNumberOfTuples(), data->getName(), true);
      T* downcastDataPtr = downcastData->getPointer(0);
      TDMSDataType::Pointer dataType = tdmsObject->dataType();
      // TODO: For now, assert that all types are double; eventually this can be made generic for all types,
      //      but we know that the data coming from the PrintRite are all stored as doubles
      assert(dataType->name() == "tdsTypeDoubleFloat");
      DoubleArrayType::Pointer tdmsData = std::dynamic_pointer_cast<DoubleArrayType>(data);
      double* tdmsDataPtr = tdmsData->getPointer(0);
      for(size_t j = 0; j < data->getNumberOfTuples(); j++)
      {
        downcastDataPtr[j] = static_cast<T>(tdmsDataPtr[j]);
      }
      downcastArrays[i] = downcastData;
    }
    return downcastArrays;
  }

private:
  std::unordered_set<std::string> m_HfChannelNames;
  std::unordered_set<std::string> m_LfChannelNames;

  std::vector<std::string> matchChannelsOfType(ChannelType type, std::unordered_map<std::string, TDMSObject::Pointer> channelObjects,
                                               const std::vector<std::string>& ignoreList = std::vector<std::string>())
  {
    std::vector<std::string> channelNames;
    std::unordered_set<std::string> localChannels = type == ChannelType::HF ? m_HfChannelNames : m_LfChannelNames;

    if(type == ChannelType::HF || type == ChannelType::LF)
    {
      for(auto&& name : localChannels)
      {
        if(std::find(std::begin(ignoreList), std::end(ignoreList), name) != std::end(ignoreList))
        {
          continue;
        }
        if(channelObjects.count(name) > 0)
        {
          channelNames.push_back(name);
        }
      }
    }
    else if(type == ChannelType::Unknown)
    {
      std::unordered_set<std::string> allChannels;
      for(auto&& it : m_HfChannelNames)
      {
        allChannels.insert(it);
      }
      for(auto&& it : m_LfChannelNames)
      {
        allChannels.insert(it);
      }

      for(auto&& pair : channelObjects)
      {
        if(std::find(std::begin(ignoreList), std::end(ignoreList), pair.first) != std::end(ignoreList) || allChannels.count(pair.first) > 0)
        {
          continue;
        }
        channelNames.push_back(pair.first);
      }
    }

    return channelNames;
  }
};

struct BoundingBox
{
  BoundingBox()
  {
    xmin = std::numeric_limits<float>::max();
    xmax = std::numeric_limits<float>::lowest();
    ymin = std::numeric_limits<float>::max();
    ymax = std::numeric_limits<float>::lowest();
  }

  bool inside(float x, float y)
  {
    bool in = false;
    if(x >= xmin && y >= ymin && x <= xmax && y <= ymax)
    {
      in = true;
    }
    return in;
  }

  void pad(float pad)
  {
    xmin -= pad;
    xmax += pad;
    ymin -= pad;
    ymax += pad;
  }

  float xmin;
  float xmax;
  float ymin;
  float ymax;
};

struct Vertex
{
  Vertex(const float& xpos, const float& ypos, const float& tol = 0.0001f)
  {
    x = xpos;
    y = ypos;
    tolerance = std::abs(tol);
  }

  float x;
  float y;
  float tolerance;
};

inline bool operator==(const Vertex& vert0, const Vertex& vert1)
{
  float tol = vert0.tolerance < vert1.tolerance ? vert0.tolerance : vert1.tolerance;
  return (SIMPLibMath::closeEnough<float>(vert0.x, vert1.x, tol) && SIMPLibMath::closeEnough<float>(vert0.y, vert1.y, tol));
}

struct Polygon
{
  Polygon() = default;

  Polygon(const std::vector<Vertex>& verts)
  {
    vertices = verts;
  }

  bool valid() const
  {
    return vertices.size() != 0;
    // return false ? vertices.empty() : true;
  }

  float signed_area() const
  {
    float area = 0.0f;
    for(auto v = std::begin(vertices); v != std::prev(std::end(vertices)); ++v)
    {
      area += ((*v).x * (*std::next(v)).y) - ((*std::next(v)).x * (*v).y);
    }
    area += (vertices.back().x * vertices.front().y) - (vertices.front().x * vertices.back().y);
    return (area * 0.5f);
  }

  std::vector<float> centroid() const
  {
    float centroidX = 0.0f;
    float centroidY = 0.0f;
    for(auto v = std::begin(vertices); v != std::prev(std::end(vertices)); ++v)
    {
      centroidX += ((*v).x + (*std::next(v)).x) * (((*v).x * (*std::next(v)).y) - ((*std::next(v)).x * (*v).y));
      centroidY += ((*v).y + (*std::next(v)).y) * (((*v).x * (*std::next(v)).y) - ((*std::next(v)).x * (*v).y));
    }
    centroidX += (vertices.back().x + vertices.front().x) * ((vertices.back().x * vertices.front().y) - (vertices.front().x * vertices.back().y));
    centroidY += (vertices.back().y + vertices.front().y) * ((vertices.back().x * vertices.front().y) - (vertices.front().x * vertices.back().y));
    centroidX /= (6.0f * signed_area());
    centroidY /= (6.0f * signed_area());
    return {centroidX, centroidY};
  }

  std::vector<float> second_moments() const
  {
    std::vector<Vertex> centeredVertices = vertices;
    std::vector<float> coa = centroid();
    std::for_each(std::begin(centeredVertices), std::end(centeredVertices), [&](Vertex& vert) {
      vert.x = vert.x - coa[0];
      vert.y = vert.y - coa[1];
    });
    float Ix = 0.0f;
    float Iy = 0.0f;
    for(auto v = std::begin(centeredVertices); v != std::prev(std::end(centeredVertices)); ++v)
    {
      Ix += (((*v).y * (*v).y) + ((*v).y * (*std::next(v)).y) + ((*std::next(v)).y * (*std::next(v)).y)) * (((*v).x * (*std::next(v)).y) - ((*std::next(v)).x * (*v).y));
      Iy += (((*v).x * (*v).x) + ((*v).x * (*std::next(v)).x) + ((*std::next(v)).x * (*std::next(v)).x)) * (((*v).x * (*std::next(v)).y) - ((*std::next(v)).x * (*v).y));
    }
    Ix += ((centeredVertices.back().y * centeredVertices.back().y) + (centeredVertices.back().y * centeredVertices.front().y) + (centeredVertices.front().y * centeredVertices.front().y)) *
          ((centeredVertices.back().x * centeredVertices.front().y) - (centeredVertices.front().x * centeredVertices.back().y));
    Iy += ((centeredVertices.back().x * centeredVertices.back().x) + (centeredVertices.back().x * centeredVertices.front().x) + (centeredVertices.front().x * centeredVertices.front().x)) *
          ((centeredVertices.back().x * centeredVertices.front().y) - (centeredVertices.front().x * centeredVertices.back().y));
    Ix /= 12.0f;
    Iy /= 12.0f;
    Ix = std::abs(Ix);
    Iy = std::abs(Iy);
    return {Ix, Iy};
  }

  std::vector<float> moment_lengths() const
  {
    std::vector<float> moments = second_moments();
    float modIx = 0.5f * std::pow((12 * 12 * moments[1] * moments[1] * moments[1] / moments[0]), 0.125f);
    float modIy = 0.5f * std::pow((12 * 12 * moments[0] * moments[0] * moments[0] / moments[1]), 0.125f);
    return {modIx, modIy};
  }

  std::vector<std::vector<float>> key_points() const
  {
    std::vector<float> coa = centroid();
    std::vector<float> ls = moment_lengths();
    std::vector<std::vector<float>> keys(2, std::vector<float>(5, 0.0f));
    keys[0][0] = coa[0];
    keys[1][0] = coa[1];
    keys[0][1] = coa[0] - ls[0];
    keys[1][1] = coa[1];
    keys[0][2] = coa[0] + ls[0];
    keys[1][2] = coa[1];
    keys[0][3] = coa[0];
    keys[1][3] = coa[1] - ls[1];
    keys[0][4] = coa[0];
    keys[1][4] = coa[1] + ls[1];
    return keys;
  }

  BoundingBox bounding_box() const
  {
    BoundingBox bbox;
    for(auto&& vert : vertices)
    {
      if(vert.x < bbox.xmin)
      {
        bbox.xmin = vert.x;
      }
      if(vert.x > bbox.xmax)
      {
        bbox.xmax = vert.x;
      }
      if(vert.y < bbox.ymin)
      {
        bbox.ymin = vert.y;
      }
      if(vert.y > bbox.ymax)
      {
        bbox.ymax = vert.y;
      }
    }
    return bbox;
  }

  std::vector<Vertex> vertices;
};

struct Polygons
{
  Polygons() = default;

  Polygons(const std::vector<Polygon>& polys)
  {
    polygons = polys;
  }

  bool exists(std::vector<Polygon>::size_type p) const
  {
    if(p >= polygons.size())
    {
      return false;
    }

    return polygons[p].valid();
  }

  std::vector<Polygon> polygons;
};

class Polynomial
{
public:
  Polynomial() = default;

  Polynomial(const std::vector<std::vector<float>>& xs, const std::vector<std::vector<float>>& ys, const uint8_t& o)
  {
    m_Inputs = xs;
    m_Output = ys;
    m_Order = o;
    m_NumDependentVariables = xs.size();
    m_NumIndependentVariables = ys.size();
    m_NumSamples = xs[0].size();
  }

  void setDependentVariables(const std::vector<std::vector<float>>& xs)
  {
    m_Inputs = xs;
    m_NumDependentVariables = xs.size();
    m_NumSamples = xs[0].size();
  }

  void setIndependentVariables(const std::vector<std::vector<float>>& ys)
  {
    m_Output = ys;
    m_NumIndependentVariables = ys.size();
  }

  void setOrder(const uint8_t& o)
  {
    m_Order = o;
    if(m_Exponents.empty())
    {
      initializeExponents();
    }
  }

  void performRegression()
  {
    // TODO: Asserting for now to enforce bicubic scaling on 2 dependent variables
    //      should eventually be made generic for other options
    //      will also need to fix set/get functions to properly update instance variables
    assert(m_NumDependentVariables == 2);
    assert(m_NumIndependentVariables == 2);
    assert(m_Order == 3);

    AMatrix A(m_NumSamples, 10);
    BVector bx(m_NumSamples);
    BVector by(m_NumSamples);

    if(m_Exponents.empty())
    {
      initializeExponents();
    }

    size_t counter = 0;
    for(auto&& exp : m_Exponents)
    {
      for(auto i = 0; i < m_NumSamples; i++)
      {
        A(i, counter) = std::pow(m_Output[0][i], static_cast<double>(exp.first)) * std::pow(m_Output[1][i], static_cast<double>(exp.second));
      }
      counter++;
    }

    for(auto i = 0; i < m_NumSamples; i++)
    {
      bx(i) = m_Inputs[0][i];
      by(i) = m_Inputs[1][i];
    }

    BVector x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bx);
    BVector y = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(by);
    m_Coefficients.push_back(x);
    m_Coefficients.push_back(y);
  }

  template <typename T>
  T transformPoint(T* pt, size_t idx)
  {
    T transformedPt = 0.0f;
    if(m_Exponents.empty())
    {
      initializeExponents();
    }
    for(auto e = 0; e < m_Exponents.size(); e++)
    {
      transformedPt += m_Coefficients[idx][e] * (std::pow(pt[0], static_cast<T>(m_Exponents[e].first)) * std::pow(pt[1], static_cast<T>(m_Exponents[e].second)));
    }
    return transformedPt;
  }

  std::vector<Eigen::VectorXd> coefficients()
  {
    return m_Coefficients;
  }

  void setCoefficients(std::vector<Eigen::VectorXd> coefficients)
  {
    // TODO: update local instance variables
    m_Coefficients = coefficients;
  }

  void setCoefficients(DoubleArrayType::Pointer coefficients)
  {
    // TODO: update local instance variables
    std::vector<BVector> localCoeffs;
    BVector x;
    BVector y;
    x.resize(10);
    y.resize(10);
    for(size_t i = 0; i < coefficients->getNumberOfTuples(); i++)
    {
      x(i) = coefficients->getComponent(i, 0);
      y(i) = coefficients->getComponent(i, 1);
    }
    localCoeffs.push_back(x);
    localCoeffs.push_back(y);
    m_Coefficients = localCoeffs;
  }

  DoubleArrayType::Pointer coefficientsAsDataArray()
  {
    std::vector<size_t> cDims(1, m_Coefficients.size());
    DoubleArrayType::Pointer coefficients = DoubleArrayType::CreateArray(m_Exponents.size(), cDims, "Spatial Scaling Coefficients", true);
    double* coefficientsPtr = coefficients->getPointer(0);
    for(auto i = 0; i < m_Exponents.size(); i++)
    {
      for(auto j = 0; j < m_Coefficients.size(); j++)
      {
        coefficientsPtr[m_Coefficients.size() * i + j] = m_Coefficients[j][i];
      }
    }
    return coefficients;
  }

  bool nullAxisCoefficients(size_t axis)
  {
    bool allNull = true;
    for(auto i = 0; i < m_Coefficients[axis].size(); i++)
    {
      if(m_Coefficients[axis](i) != 0)
      {
        allNull = false;
        break;
      }
    }
    return allNull;
  }

  bool nullCoefficients()
  {
    bool allNull = true;
    for(auto i = 0; i < m_Coefficients.size(); i++)
    {
      allNull = nullAxisCoefficients(i);
      if(!allNull)
      {
        break;
      }
    }
    return allNull;
  }

private:
  void initializeExponents()
  {
    for(uint8_t i = 0; i <= m_Order; i++)
    {
      for(uint8_t j = 0; j <= m_Order; j++)
      {
        if((i + j) <= m_Order)
        {
          m_Exponents.emplace_back(std::make_pair(j, i));
        }
      }
    }
  }

  using AMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  using BVector = Eigen::VectorXd;

  std::vector<std::vector<float>> m_Inputs;
  std::vector<std::vector<float>> m_Output;
  uint8_t m_Order;
  std::vector<std::vector<float>>::size_type m_NumDependentVariables;
  std::vector<std::vector<float>>::size_type m_NumIndependentVariables;
  std::vector<float>::size_type m_NumSamples;
  std::vector<std::pair<size_t, size_t>> m_Exponents;
  std::vector<BVector> m_Coefficients;
};
}; // namespace PrintRiteHelpers
