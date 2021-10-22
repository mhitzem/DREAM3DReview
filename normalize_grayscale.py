import numpy as np
import os

from contextlib import ExitStack
from pathlib import Path
from typing import Any, AnyStr, Dict, Generator, List, Match, Pattern
from enum import IntEnum
from typing import List, Tuple, Union

from dream3d.Filter import Filter, FilterDelegatePy
import dream3d.simpl as simpl
import dream3d.simplpy as simplpy

class NormalizeGrayscale(Filter):
  def __init__(self) -> None:
    self.data_array: DataArrayPath = simpl.DataArrayPath('', '', '')
    self.arr_name: str = ''
    self.choice: bool = False
    self.axes: int = 0
    self.new_arr: DataArrayPath = simpl.DataArrayPath('', '', '')

  def _set_data_array(self, value: simpl.DataArrayPath) -> None:
    self.data_array = value

  def _get_data_array(self) -> simpl.DataArrayPath:
    return self.data_array

  def _set_axes(self, value: int) -> None:
    self.axes = value

  def _get_axes(self) -> int:
    return self.axes

  def _set_choice(self, value: bool) -> None:
    self.choice = value

  def _get_choice(self) -> bool:
    return self.choice

  def _set_arr_name(self, value: str) -> None:
    self.arr_name = value

  def _get_arr_name(self) -> str:
    return self.arr_name

  def _set_new_arr(self, value: simpl.DataArrayPath) -> None:
    self.new_arr = value

  def _get_new_arr(self) -> simpl.DataArrayPath:
    return self.new_arr

  @staticmethod
  def name() -> str:
    return 'Grayscale Normalizer'

  @staticmethod
  def uuid() -> str:
    return '{b98fa052-4c74-4d25-96ce-5b95074fcecf}'

  @staticmethod
  def group_name() -> str:
    return 'Example'

  @staticmethod
  def sub_group_name() -> str:
    return 'Sub Example'

  @staticmethod
  def human_label() -> str:
    return 'Normalize Grayscale Values'

  @staticmethod
  def version() -> str:
    return '1.0.0'

  @staticmethod
  def compiled_lib_name() -> str:
    return 'Python'

  def setup_parameters(self) -> List[simpl.FilterParameter]:
    req = simpl.DataArraySelectionFilterParameter.RequirementType([simpl.IGeometry.Type.Image], [], [], [])
    return [
      simpl.DataArraySelectionFilterParameter('Select Data Array', 'data_array', self.data_array, simpl.FilterParameter.Category.RequiredArray, self._set_data_array, self._get_data_array, req, -1),
      simpl.ChoiceFilterParameter('Slice Axis', 'axes', self.axes, simpl.FilterParameter.Category.Parameter, self._set_axes,
                                 self._get_axes, ["x", "y", "z"], False, -1),
      simpl.LinkedBooleanFilterParameter('Create New Array', 'create_new_array', self.choice, simpl.FilterParameter.Category.Parameter, self._set_choice, self._get_choice, ['new_array'], -1),
      simpl.DataArrayCreationFilterParameter('New Array', 'new_array', self.new_arr, simpl.FilterParameter.Category.CreatedArray, self._set_new_arr, self._get_new_arr, simpl.DataArrayCreationFilterParameter.RequirementType(), -1)
    ]

  def data_check(self, dca: simpl.DataContainerArray, status_delegate: Union[simpl.FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    dc = dca.getDataContainer(self.data_array)
    am = dca.getAttributeMatrix(self.data_array)
    selected_data_array = am.getAttributeArray(self.data_array.DataArrayName)
    #selected_data_array.getNumComponents != 1 -> could be req

    if not dca.doesAttributeMatrixExist(self.data_array):
      return (-5550, 'One data array must be selected')

    if type(dc.Geometry.getDimensions()) != simpl.SizeVec3:
      return (-301, 'Not a 3D array')

    if not dc.Geometry:
      return (-302, 'DataContainer has no geometry')

    if dc.Geometry.getGeometryType() != simpl.IGeometry.Type.Image:
      return (-303, 'Wrong geometry type')

    if self.choice:
      udims = dc.Geometry.getDimensions()
      shape = [np.int64(udims[2]), np.int64(udims[1]), np.int64(udims[0])]
      arr_dtype = am.getAttributeArray(self.data_array).dtype

      new_data_array = np.zeros(shape, dtype = arr_dtype)

    # Create a simpl.DataArray object to hold the vertices by reference
      #switch statement for dtype
      new_simpl_array = simpl.UInt8Array(new_data_array, self.new_arr.DataArrayName)
      am.addOrReplaceAttributeArray(new_simpl_array)

    return (0, 'Success')

  def _execute_impl(self, dca: simpl.DataContainerArray, status_delegate: Union[simpl.FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    dc = dca.getDataContainer(self.data_array)
    udims = dc.Geometry.getDimensions()
    shape = [np.int64(udims[2]), np.int64(udims[1]), np.int64(udims[0])]

    da = dca.getAttributeMatrix(self.data_array).getAttributeArray(self.data_array)
    data = da.npview()
    data = np.reshape(data, shape)

    if self.choice:
      #simpl_data_array
      da_new = dca.getAttributeMatrix(self.new_arr).getAttributeArray(self.new_arr)
      #npv_data_array
      corrected_data = da_new.npview()
      corrected_data = np.reshape(corrected_data, shape)

    arr = []

    slice_axis = self.axes
    shape = data.shape

    arr_avg = np.average(data)

    #slices are in x-axis
    if slice_axis == 0:
      layer_avg = np.average(data, axis=(1,2))
      correction_factor = layer_avg/arr_avg
      corrected_data = data / correction_factor[:, None, None]

    #slices are in y-axis
    elif slice_axis == 1:
      layer_avg = np.average(data, axis=(0,2))
      correction_factor = layer_avg / arr_avg
      corrected_data = data/correction_factor[None, :, None]

    #slices are in z-axis
    elif slice_axis == 2:
      layer_avg = np.average(data, axis=(0,1))
      correction_factor = layer_avg / arr_avg
      corrected_data = data / correction_factor[None, None, :]

    if not self.choice:
      data = corrected_data

    return (0, 'Success')


filters = [NormalizeGrayscale]
