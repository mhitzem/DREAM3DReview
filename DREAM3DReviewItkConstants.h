/* ============================================================================
 * Copyright (c) 2021-2021 BlueQuartz Software, LLC
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
 *    United States Air Force Prime Contract FA8650-10-D-5210
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#if defined(SIMPL_USE_ITK)

#include "itkConfigure.h"

#if defined(ITK_VERSION_MAJOR) && ITK_VERSION_MAJOR == 4
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wself-assign-field"
#endif
#endif

#include "itkImage.h"
#include "itkRGBPixel.h"
#endif

/**
 * @brief This namespace is used to define some Constants for the plugin itself.
 */
namespace DREAM3DReviewConstants
{
#if defined(SIMPL_USE_ITK)

// multicomponent pixels
using RGBUInt8PixelType = itk::RGBPixel<uint8_t>; // ipf color etc
// typedef itk::RGBAPixel <float> RGBAFloatPixelType; //may be able to handle quats with this?

// define pixels for dream3d variable types
using Int8PixelType = int8_t;
using UInt8PixelType = uint8_t;
using Int16PixelType = int16_t;
using UInt16PixelType = uint16_t;
using Int32PixelType = int32_t;
using UInt32PixelType = uint32_t;
using Int64PixelType = int64_t;
using UInt64PixelType = uint64_t;
using FloatPixelType = float;
using DoublePixelType = double;

// define default pixel type
#if Anisotropy_BitDepth == 8
typedef UInt8PixelType DefaultPixelType;
typedef DataArray<DefaultPixelType> DefaultArrayType;
#elif Anisotropy_BitDepth == 16
typedef UInt16PixelType DefaultPixelType;
typedef UInt16ArrayType DefaultArrayType;
#elif Anisotropy_BitDepth == 32
typedef FloatPixelType DefaultPixelType;
typedef FloatArrayType DefaultArrayType;
#else
using DefaultPixelType = UInt8PixelType;
using DefaultArrayType = UInt8ArrayType;
#endif

// define dimensionality
const unsigned int SliceDimension = 2;
const unsigned int ImageDimension = 3;

// slice directions
const unsigned int XSlice = 0;
const unsigned int YSlice = 1;
const unsigned int ZSlice = 2;

// define image types
using DefaultImageType = itk::Image<DefaultPixelType, ImageDimension>;
using Int8ImageType = itk::Image<Int8PixelType, ImageDimension>;
using UInt8ImageType = itk::Image<UInt8PixelType, ImageDimension>;
using Int16ImageType = itk::Image<Int16PixelType, ImageDimension>;
using UInt16ImageType = itk::Image<UInt16PixelType, ImageDimension>;
using Int32ImageType = itk::Image<Int32PixelType, ImageDimension>;
using UInt32ImageType = itk::Image<UInt32PixelType, ImageDimension>;
using Int64ImageType = itk::Image<Int64PixelType, ImageDimension>;
using UInt64ImageType = itk::Image<UInt64PixelType, ImageDimension>;
using FloatImageType = itk::Image<FloatPixelType, ImageDimension>;
using DoubleImageType = itk::Image<DoublePixelType, ImageDimension>;

using DefaultSliceType = itk::Image<DefaultPixelType, SliceDimension>;
using Int8SliceType = itk::Image<Int8PixelType, SliceDimension>;
using UInt8SliceType = itk::Image<UInt8PixelType, SliceDimension>;
using Int16SliceType = itk::Image<Int16PixelType, SliceDimension>;
using UInt16SliceType = itk::Image<UInt16PixelType, SliceDimension>;
using Int32SliceType = itk::Image<Int32PixelType, SliceDimension>;
using UInt32SliceType = itk::Image<UInt32PixelType, SliceDimension>;
using Int64SliceType = itk::Image<Int64PixelType, SliceDimension>;
using UInt64SliceType = itk::Image<UInt64PixelType, SliceDimension>;
using FloatSliceType = itk::Image<FloatPixelType, SliceDimension>;
using DoubleSliceType = itk::Image<DoublePixelType, SliceDimension>;
#endif
} // namespace DREAM3DReviewConstants
