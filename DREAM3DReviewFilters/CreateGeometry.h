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
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef _creategeometry_h_
#define _creategeometry_h_

#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The CreateGeometry class. See [Filter documentation](@ref creategeometry) for details.
 */
class CreateGeometry : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(CreateGeometry)
  SIMPL_STATIC_NEW_MACRO(CreateGeometry)
  SIMPL_TYPE_MACRO_SUPER(CreateGeometry, AbstractFilter)

  virtual ~CreateGeometry();

  SIMPL_FILTER_PARAMETER(int, GeometryType)
  Q_PROPERTY(int GeometryType READ getGeometryType WRITE setGeometryType)

  SIMPL_FILTER_PARAMETER(QString, DataContainerName)
  Q_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedVertexListArrayPath0)
  Q_PROPERTY(DataArrayPath SharedVertexListArrayPath0 READ getSharedVertexListArrayPath0 WRITE setSharedVertexListArrayPath0)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedVertexListArrayPath1)
  Q_PROPERTY(DataArrayPath SharedVertexListArrayPath1 READ getSharedVertexListArrayPath1 WRITE setSharedVertexListArrayPath1)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedVertexListArrayPath2)
  Q_PROPERTY(DataArrayPath SharedVertexListArrayPath2 READ getSharedVertexListArrayPath2 WRITE setSharedVertexListArrayPath2)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedVertexListArrayPath3)
  Q_PROPERTY(DataArrayPath SharedVertexListArrayPath3 READ getSharedVertexListArrayPath3 WRITE setSharedVertexListArrayPath3)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedVertexListArrayPath4)
  Q_PROPERTY(DataArrayPath SharedVertexListArrayPath4 READ getSharedVertexListArrayPath4 WRITE setSharedVertexListArrayPath4)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedEdgeListArrayPath)
  Q_PROPERTY(DataArrayPath SharedEdgeListArrayPath READ getSharedEdgeListArrayPath WRITE setSharedEdgeListArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedTriListArrayPath)
  Q_PROPERTY(DataArrayPath SharedTriListArrayPath READ getSharedTriListArrayPath WRITE setSharedTriListArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedQuadListArrayPath)
  Q_PROPERTY(DataArrayPath SharedQuadListArrayPath READ getSharedQuadListArrayPath WRITE setSharedQuadListArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, SharedTetListArrayPath)
  Q_PROPERTY(DataArrayPath SharedTetListArrayPath READ getSharedTetListArrayPath WRITE setSharedTetListArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, XBoundsArrayPath)
  Q_PROPERTY(DataArrayPath XBoundsArrayPath READ getXBoundsArrayPath WRITE setXBoundsArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, YBoundsArrayPath)
  Q_PROPERTY(DataArrayPath YBoundsArrayPath READ getYBoundsArrayPath WRITE setYBoundsArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, ZBoundsArrayPath)
  Q_PROPERTY(DataArrayPath ZBoundsArrayPath READ getZBoundsArrayPath WRITE setZBoundsArrayPath)

  SIMPL_FILTER_PARAMETER(IntVec3_t, Dimensions)
  Q_PROPERTY(IntVec3_t Dimensions READ getDimensions WRITE setDimensions)

  SIMPL_FILTER_PARAMETER(FloatVec3_t, Origin)
  Q_PROPERTY(FloatVec3_t Origin READ getOrigin WRITE setOrigin)

  SIMPL_FILTER_PARAMETER(FloatVec3_t, Resolution)
  Q_PROPERTY(FloatVec3_t Resolution READ getResolution WRITE setResolution)

  SIMPL_FILTER_PARAMETER(QString, ImageCellAttributeMatrixName)
  Q_PROPERTY(QString ImageCellAttributeMatrixName READ getImageCellAttributeMatrixName WRITE setImageCellAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, RectGridCellAttributeMatrixName)
  Q_PROPERTY(QString RectGridCellAttributeMatrixName READ getRectGridCellAttributeMatrixName WRITE setRectGridCellAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName0)
  Q_PROPERTY(QString VertexAttributeMatrixName0 READ getVertexAttributeMatrixName0 WRITE setVertexAttributeMatrixName0)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName1)
  Q_PROPERTY(QString VertexAttributeMatrixName1 READ getVertexAttributeMatrixName1 WRITE setVertexAttributeMatrixName1)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName2)
  Q_PROPERTY(QString VertexAttributeMatrixName2 READ getVertexAttributeMatrixName2 WRITE setVertexAttributeMatrixName2)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName3)
  Q_PROPERTY(QString VertexAttributeMatrixName3 READ getVertexAttributeMatrixName3 WRITE setVertexAttributeMatrixName3)

  SIMPL_FILTER_PARAMETER(QString, VertexAttributeMatrixName4)
  Q_PROPERTY(QString VertexAttributeMatrixName4 READ getVertexAttributeMatrixName4 WRITE setVertexAttributeMatrixName4)

  SIMPL_FILTER_PARAMETER(QString, EdgeAttributeMatrixName)
  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, FaceAttributeMatrixName0)
  Q_PROPERTY(QString FaceAttributeMatrixName0 READ getFaceAttributeMatrixName0 WRITE setFaceAttributeMatrixName0)

  SIMPL_FILTER_PARAMETER(QString, FaceAttributeMatrixName1)
  Q_PROPERTY(QString FaceAttributeMatrixName1 READ getFaceAttributeMatrixName1 WRITE setFaceAttributeMatrixName1)

  SIMPL_FILTER_PARAMETER(QString, TetCellAttributeMatrixName)
  Q_PROPERTY(QString TetCellAttributeMatrixName READ getTetCellAttributeMatrixName WRITE setTetCellAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(bool, TreatWarningsAsErrors)
  Q_PROPERTY(bool TreatWarningsAsErrors READ getTreatWarningsAsErrors WRITE setTreatWarningsAsErrors)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  virtual const QString getCompiledLibraryName();

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
  */
  virtual const QString getBrandingString();

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  virtual const QString getFilterVersion();

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  virtual const QString getGroupName();

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  virtual const QString getSubGroupName();

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  virtual const QString getHumanLabel();

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  virtual void setupFilterParameters();

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  virtual void readFilterParameters(AbstractFilterParametersReader* reader, int index);

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  virtual void execute();

  /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  virtual void preflight();

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  CreateGeometry();
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(float, XBounds)
  DEFINE_DATAARRAY_VARIABLE(float, YBounds)
  DEFINE_DATAARRAY_VARIABLE(float, ZBounds)
  DEFINE_DATAARRAY_VARIABLE(int64_t, Edges)
  DEFINE_DATAARRAY_VARIABLE(int64_t, Tris)
  DEFINE_DATAARRAY_VARIABLE(int64_t, Quads)
  DEFINE_DATAARRAY_VARIABLE(int64_t, Tets)

  int64_t m_NumVerts;

  CreateGeometry(const CreateGeometry&); // Copy Constructor Not Implemented
  void operator=(const CreateGeometry&); // Operator '=' Not Implemented
};

#endif /* _CreateGeometry_H_ */
