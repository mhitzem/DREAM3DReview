/* ============================================================================
 * Copyright (c) 2009-2016 BlueQuartz Software, LLC
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
 * contributors may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The code contained herein was partially funded by the following contracts:
 *    United States Air Force Prime Contract FA8650-07-D-5800
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *    United States Prime Contract Navy N00173-07-C-2068
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "NormalizeArrays.h"

#include <cassert>
#include <cmath>

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

#include <QtCore/QTextStream>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/IDataArray.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedChoicesFilterParameter.h"
#include "SIMPLib/FilterParameters/MultiDataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

/* Create Enumerations to allow the created Attribute Arrays to take part in renaming */
enum createdPathID : RenameDataPath::DataID_t
{
  DataArrayID30 = 30,
  DataArrayID31 = 31,
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NormalizeArrays::NormalizeArrays() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
NormalizeArrays::~NormalizeArrays() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::setupFilterParameters()
{
  FilterParameterVectorType parameters;
  {
    LinkedChoicesFilterParameter::Pointer parameter = LinkedChoicesFilterParameter::New();
    parameter->setHumanLabel("Operation Type");
    parameter->setPropertyName("NormalizeType");
    parameter->setSetterCallback(SIMPL_BIND_SETTER(NormalizeArrays, this, NormalizeType));
    parameter->setGetterCallback(SIMPL_BIND_GETTER(NormalizeArrays, this, NormalizeType));
    std::vector<QString> choices;
    choices.push_back("Rescale to Range");
    choices.push_back("Standardize");
    parameter->setChoices(choices);
    std::vector<QString> linkedProps = {"RangeMin", "RangeMax"};
    parameter->setLinkedProperties(linkedProps);
    parameter->setEditable(false);
    parameter->setCategory(FilterParameter::Category::Parameter);
    parameters.push_back(parameter);
  }
  parameters.push_back(SIMPL_NEW_STRING_FP("Postfix", Postfix, FilterParameter::Category::Parameter, NormalizeArrays));
  std::vector<QString> linkedProps = {"MaskArrayPath", "DefaultValue"};
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Category::Parameter, NormalizeArrays, linkedProps));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Default Masked Value", DefaultValue, FilterParameter::Category::Parameter, NormalizeArrays));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Range Minimum", RangeMin, FilterParameter::Category::Parameter, NormalizeArrays, 0));
  parameters.push_back(SIMPL_NEW_DOUBLE_FP("Range Maximum", RangeMax, FilterParameter::Category::Parameter, NormalizeArrays, 0));
  MultiDataArraySelectionFilterParameter::RequirementType req =
      MultiDataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_MDA_SELECTION_FP("Attribute Arrays to Normalize", SelectedDataArrayPaths, FilterParameter::Category::RequiredArray, NormalizeArrays, req));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::Category::RequiredArray, NormalizeArrays, dasReq));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::initialize()
{
  m_SelectedWeakPtrVector.clear();
  m_NormalizedArraysPtrVector.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::dataCheck()
{
  clearErrorCode();
  clearWarningCode();
  initialize();

  if(getNormalizeType() != 0 && getNormalizeType() != 1)
  {
    QString ss = QObject::tr("Invalid selection for operation type");
    setErrorCondition(-701, ss);
    return;
  }

  if(getSelectedDataArrayPaths().empty())
  {
    QString ss = QObject::tr("At least one Attribute Array must be selected");
    setErrorCondition(-11001, ss);
    return;
  }

  if(getPostfix().isEmpty())
  {
    QString ss = QObject::tr("A postfix for the normalized Attribute Arrays must be entered");
    setErrorCondition(-11001, ss);
  }

  std::vector<DataArrayPath> paths = getSelectedDataArrayPaths();

  if(!DataArrayPath::ValidateVector(paths))
  {
    QString ss = QObject::tr("There are Attribute Arrays selected that are not contained in the same Attribute Matrix; all selected Attribute Arrays must belong to the same Attribute Matrix");
    setErrorCondition(-11002, ss);
  }

  if(getErrorCode() < 0)
  {
    return;
  }

  std::vector<size_t> cDims(1, 1);
  QString dcName = paths[0].getDataContainerName();
  QString amName = paths[0].getAttributeMatrixName();

  for(auto&& path : paths)
  {
    IDataArray::WeakPointer ptr = getDataContainerArray()->getPrereqIDataArrayFromPath(this, path);
    if(ptr.lock())
    {
      m_SelectedWeakPtrVector.push_back(ptr);
      int32_t numComps = ptr.lock()->getNumberOfComponents();
      if(numComps != 1)
      {
        QString ss = QObject::tr("All Attribute Arrays must be scalar arrays, but %1 has %2 total components").arg(ptr.lock()->getName()).arg(numComps);
        setErrorCondition(-11003, ss);
      }
      else
      {
        QString arrayName = path.getDataArrayName() + getPostfix();
        DataArrayPath tempPath(dcName, amName, arrayName);
        double defaultValue = m_UseMask ? getDefaultValue() : 0.0;

        DoubleArrayType::WeakPointer ptr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<double>>(this, tempPath, defaultValue, cDims, "", DataArrayID31);
        if(getErrorCode() >= 0)
        {
          m_NormalizedArraysPtrVector.push_back(ptr.lock());
        }
      }
    }
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCode() >= 0)
    {
      paths.push_back(getMaskArrayPath());
    }
  }

  getDataContainerArray()->validateNumberOfTuples(this, paths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void copyDataArrays(IDataArray::Pointer dataPtr, std::vector<double>& copy, bool useMask, bool* mask)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  T* dPtr = inDataPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        copy.push_back(dPtr[i]);
      }
    }
    else
    {
      copy.push_back(dPtr[i]);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
class NormalizeArraysImpl
{
public:
  NormalizeArraysImpl(std::vector<std::vector<double>>& arrays, int32_t normalizeType, double rangeMin, double rangeMax)
  : m_Arrays(arrays)
  , m_NormalizeType(normalizeType)
  , m_RangeMin(rangeMin)
  , m_RangeMax(rangeMax)
  {
  }

  virtual ~NormalizeArraysImpl() = default;

  void compute(size_t start, size_t end) const
  {
    for(size_t i = start; i < end; i++)
    {
      if(m_NormalizeType == 0)
      {
        double min = *std::min_element(std::begin(m_Arrays[i]), std::end(m_Arrays[i]));
        double max = *std::max_element(std::begin(m_Arrays[i]), std::end(m_Arrays[i]));
        // size_t count = m_Arrays.size();
        for(auto& val : m_Arrays[i])
        {
          val = m_RangeMin + (((val - min) * (m_RangeMax - m_RangeMin)) / (max - min));
        }
      }
      else if(m_NormalizeType == 1)
      {
        double sum = std::accumulate(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), 0.0);
        double mean = sum / m_Arrays[i].size();
        std::vector<double> difference(m_Arrays[i].size());
        std::transform(std::begin(m_Arrays[i]), std::end(m_Arrays[i]), std::begin(difference), [&](double val) { return val - mean; });
        double squaredSum = std::inner_product(std::begin(difference), std::end(difference), std::begin(difference), 0.0f);
        double stdDev = std::sqrt(squaredSum / m_Arrays[i].size());
        for(auto& val : m_Arrays[i])
        {
          val = (val - mean) / stdDev;
        }
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    compute(r.begin(), r.end());
  }
#endif

private:
  std::vector<std::vector<double>>& m_Arrays;
  int32_t m_NormalizeType;
  double m_RangeMin;
  double m_RangeMax;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void NormalizeArrays::execute()
{
  dataCheck();
  if(getErrorCode() < 0)
  {
    return;
  }

  if(m_SelectedDataArrayPaths.size() != m_SelectedWeakPtrVector.size())
  {
    QString ss = QObject::tr("The number of selected Attribute Arrays does not equal the number of internal weak pointers");
    setErrorCondition(-11008, ss);
    return;
  }

  size_t numTuples = m_SelectedWeakPtrVector[0].lock()->getNumberOfTuples();

  std::vector<std::vector<double>> arrays(m_SelectedWeakPtrVector.size());
  for(auto&& vec : arrays)
  {
    vec.reserve(numTuples);
  }

  for(std::vector<double>::size_type i = 0; i < arrays.size(); i++)
  {
    EXECUTE_FUNCTION_TEMPLATE(this, copyDataArrays, m_SelectedWeakPtrVector[i].lock(), m_SelectedWeakPtrVector[i].lock(), arrays[i], m_UseMask, m_Mask)
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  if(true)
  {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, arrays.size()), NormalizeArraysImpl(arrays, m_NormalizeType, m_RangeMin, m_RangeMax), tbb::auto_partitioner());
  }
  else
#endif
  {
    NormalizeArraysImpl serial(arrays, m_NormalizeType, m_RangeMin, m_RangeMax);
    serial.compute(0, arrays.size());
  }

  size_t index = 0;

  assert(m_NormalizedArraysPtrVector.size() == arrays.size());

  for(size_t i = 0; i < numTuples; i++)
  {
    for(std::vector<double>::size_type j = 0; j < arrays.size(); j++)
    {
      if(m_UseMask)
      {
        if(m_Mask[i])
        {
          m_NormalizedArraysPtrVector[j]->setValue(i, arrays[j][index]);
          index++;
        }
      }
      else
      {
        m_NormalizedArraysPtrVector[j]->setValue(i, arrays[j][i]);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer NormalizeArrays::newFilterInstance(bool copyFilterParameters) const
{
  NormalizeArrays::Pointer filter = NormalizeArrays::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QUuid NormalizeArrays::getUuid() const
{
  return QUuid("{8c584519-15c3-5010-b5ed-a2ac626591a1}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString NormalizeArrays::getHumanLabel() const
{
  return "Normalize Attribute Arrays";
}

// -----------------------------------------------------------------------------
NormalizeArrays::Pointer NormalizeArrays::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
std::shared_ptr<NormalizeArrays> NormalizeArrays::New()
{
  struct make_shared_enabler : public NormalizeArrays
  {
  };
  std::shared_ptr<make_shared_enabler> val = std::make_shared<make_shared_enabler>();
  val->setupFilterParameters();
  return val;
}

// -----------------------------------------------------------------------------
QString NormalizeArrays::getNameOfClass() const
{
  return QString("NormalizeArrays");
}

// -----------------------------------------------------------------------------
QString NormalizeArrays::ClassName()
{
  return QString("NormalizeArrays");
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setSelectedDataArrayPaths(const std::vector<DataArrayPath>& value)
{
  m_SelectedDataArrayPaths = value;
}

// -----------------------------------------------------------------------------
std::vector<DataArrayPath> NormalizeArrays::getSelectedDataArrayPaths() const
{
  return m_SelectedDataArrayPaths;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setNormalizeType(int value)
{
  m_NormalizeType = value;
}

// -----------------------------------------------------------------------------
int NormalizeArrays::getNormalizeType() const
{
  return m_NormalizeType;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setRangeMin(double value)
{
  m_RangeMin = value;
}

// -----------------------------------------------------------------------------
double NormalizeArrays::getRangeMin() const
{
  return m_RangeMin;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setRangeMax(double value)
{
  m_RangeMax = value;
}

// -----------------------------------------------------------------------------
double NormalizeArrays::getRangeMax() const
{
  return m_RangeMax;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setPostfix(const QString& value)
{
  m_Postfix = value;
}

// -----------------------------------------------------------------------------
QString NormalizeArrays::getPostfix() const
{
  return m_Postfix;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setUseMask(bool value)
{
  m_UseMask = value;
}

// -----------------------------------------------------------------------------
bool NormalizeArrays::getUseMask() const
{
  return m_UseMask;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setMaskArrayPath(const DataArrayPath& value)
{
  m_MaskArrayPath = value;
}

// -----------------------------------------------------------------------------
DataArrayPath NormalizeArrays::getMaskArrayPath() const
{
  return m_MaskArrayPath;
}

// -----------------------------------------------------------------------------
void NormalizeArrays::setDefaultValue(double value)
{
  m_DefaultValue = value;
}

// -----------------------------------------------------------------------------
double NormalizeArrays::getDefaultValue() const
{
  return m_DefaultValue;
}
