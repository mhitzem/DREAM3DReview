"""
AttributeMatrixCreationFilterParameter
AttributeMatrixSelectionFilterParameter
AxisAngleFilterParameter
+ BooleanFilterParameter
CalculatorFilterParameter
ChoiceFilterParameter
ComparisonSelectionAdvancedFilterParameter
ComparisonSelectionFilterParameter
ConstrainedDoubleFilterParameter
ConstrainedIntFilterParameter
DataArrayCreationFilterParameter
+ DataArraySelectionFilterParameter
DataContainerArrayProxyFilterParameter
DataContainerCreationFilterParameter
DataContainerGridSelectionFilterParameter
DataContainerReaderFilterParameter
DataContainerSelectionFilterParameter
DoubleFilterParameter
DynamicChoiceFilterParameter
DynamicTableFilterParameter
FileListInfoFilterParameter
FloatFilterParameter
FloatVec2FilterParameter
FloatVec3FilterParameter
FourthOrderPolynomialFilterParameter
GenerateColorTableFilterParameter
ImportHDF5DatasetFilterParameter
+ InputFileFilterParameter
InputPathFilterParameter
* IntFilterParameter
IntVec2FilterParameter
IntVec3FilterParameter
LinkedBooleanFilterParameter
LinkedChoicesFilterParameter
LinkedDataContainerSelectionFilterParameter
MontageSelectionFilterParameter
MontageStructureSelectionFilterParameter
MultiAttributeMatrixSelectionFilterParameter
MultiDataContainerSelectionFilterParameter
MultiDataArraySelectionFilterParameter
NumericTypeFilterParameter
OutputFileFilterParameter
OutputPathFilterParameter
ParagraphFilterParameter
PreflightUpdatedValueFilterParameter
RangeFilterParameter
ReadASCIIDataFilterParameter
ScalarTypeFilterParameter
SecondOrderPolynomialFilterParameter
SeparatorFilterParameter
ShapeTypeSelectionFilterParameter
+ StringFilterParameter
ThirdOrderPolynomialFilterParameter
UInt64FilterParameter
UnknownFilterParameter
MultiInputFileFilterParameter"""


from typing import List, Tuple, Union
from dream3d.Filter import Filter, FilterDelegatePy
from dream3d.simpl import BooleanFilterParameter, DataContainerArray, StringFilterParameter, InputFileFilterParameter, FilterDelegateCpp, FilterParameter, IntFilterParameter, DataArraySelectionFilterParameter, DataArrayPath
from dream3d.simpl import InputPathFilterParameter, FloatFilterParameter, StackFileListInfo, FileListInfoFilterParameter

class D3DReviewTestFilter(Filter):
  def __init__(self) -> None:
    self.int_param: int = 5
    self.dap_param: DataArrayPath = DataArrayPath('', '', '')
    self.str_param: str = "Something"
    self.bool_param:bool = False
    self.input_file_param:str = "/No/Path/anywhere.txt"
    self.input_path_param:str = "/Some/Folder/Somewhere"
    self.float_param:float = 33.44
    # The arguments to StackFileListInfo are the following:
    #param paddingDigits The number of padding digits use to generate the file names
    #param ordering How to order the files: Low to High = 0; High to Low = 1
    #param startIndex The starting index
    #param endIndex The ending index (inclusive)
    #param incrementIndex How much to increment index each time through the loop
    #param inputPath The initial input file path (String)
    #param filePrefix The File prefix (String)
    #param fileSuffix The File Suffix (String)
    #param fileExtension The file extnsion without the '.' chacter

    self.stack_info_param:StackFileListInfo = StackFileListInfo(0,0,0,0,1, 'InputPath', 'File Prefix_', '_Suffix', 'tdms')

  def _set_float_param(self, value: float) -> None:
      self.float_param = value
  def _get_float_param(self) -> float:
      return self.float_param

  def _set_int(self, value: int) -> None:
    self.int_param = value

  def _get_int(self) -> int:
    return self.int_param

  def _set_dap(self, value: DataArrayPath) -> None:
    self.dap_param = value

  def _get_dap(self) -> DataArrayPath:
    return self.dap_param

  def _get_str(self) -> str:
      return self.str_param

  def _set_str(self, value: str) -> None:
      self.str_param = value

  def _get_bool(self) -> bool:
      return self.bool_param

  def _set_bool(self, value ) -> None:
      self.bool_param = value

  def _get_input_file(self) -> str:
      return self.input_file_param

  def _set_input_file(self, value: str ) -> None:
      self.input_file_param = value
  
  def _get_input_path_param(self) -> str:
      return self.input_path_param  

  def _set_input_path_param(self, value: str) -> None:
      self.input_path_param = value

  def _get_stack_info_param(self) -> StackFileListInfo:
      return self.stack_info_param

  def _set_stack_info_param(self, value) -> None:
      self.stack_info_param = value
  
  @staticmethod
  def name() -> str:
    return 'D3DReviewTestFilter'

  @staticmethod
  def uuid() -> str:
    return '{6128bd81-8c49-54c4-b4c7-06fab6ae456a}'

  @staticmethod
  def group_name() -> str:
    return 'Test'

  @staticmethod
  def sub_group_name() -> str:
    return 'Test'

  @staticmethod
  def human_label() -> str:
    return '_Test Filter [Python]'

  @staticmethod
  def version() -> str:
    return '0.1.0'

  @staticmethod
  def compiled_lib_name() -> str:
    return 'DREAM3DReview [Python]'

  def setup_parameters(self) -> List[FilterParameter]:
    req = DataArraySelectionFilterParameter.RequirementType()
    return [
        IntFilterParameter('Integer', 'int_param', self.int_param, FilterParameter.Category.Parameter, self._set_int, self._get_int, -1),
        DataArraySelectionFilterParameter('Data Array Path Selection', 'dap_param', self.dap_param, FilterParameter.Category.RequiredArray, self._set_dap, self._get_dap, req, -1),
        StringFilterParameter('String', 'str_param', self.str_param, FilterParameter.Category.Parameter, self._set_str, self._get_str, -1),
        BooleanFilterParameter('Boolean', 'bool_param', self.bool_param, FilterParameter.Category.Parameter, self._set_bool, self._get_bool, -1),
        InputFileFilterParameter('Input File', 'input_file_param', self.input_file_param, 
                                    FilterParameter.Category.Parameter, self._set_input_file, self._get_input_file,'*.ang', 'EDAX Ang', -1),
        InputPathFilterParameter('Input Directory', 'input_path_param', self.input_path_param, 
                                    FilterParameter.Category.Parameter, self._set_input_path_param, self._get_input_path_param, -1),
        FileListInfoFilterParameter('List of TDMS Files', 'stack_info_param', self.stack_info_param,
                                    FilterParameter.Category.Parameter, self._set_stack_info_param, self._get_stack_info_param)
    ]

  def data_check(self, dca: DataContainerArray, delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    am = dca.getAttributeMatrix(self.dap_param)

    if am is None:
      return (-1, 'AttributeMatrix is None')

    da = am.getAttributeArray(self.dap_param)
    if da is None:
      return (-2, 'DataArray is None')

    delegate.notifyStatusMessage('data_check finished!')

    return (0, 'Success')

  def _execute_impl(self, dca: DataContainerArray, delegate: Union[FilterDelegateCpp, FilterDelegatePy] = FilterDelegatePy()) -> Tuple[int, str]:
    delegate.notifyStatusMessage(f'int_param = {self.int_param}')

    da = dca.getAttributeMatrix(self.dap_param).getAttributeArray(self.dap_param)

    data = da.npview()
    delegate.notifyStatusMessage(f'before = {data}')
    data += self.int_param
    delegate.notifyStatusMessage(f'after = {data}')

    delegate.notifyStatusMessage('execute finished!')
    return (0, 'Success')

filters = [D3DReviewTestFilter]