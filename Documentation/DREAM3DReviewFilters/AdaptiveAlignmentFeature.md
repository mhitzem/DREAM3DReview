# Adaptive Alignment (Feature)  #


## Group (Subgroup) ##

Unsupported (Anisotropy)

## Description ##

This **Filter** attempts to align 'sections' of the sample perpendicular to the Z-direction by determining the position that results in the most overlap of previously defined "regions".  The "regions" are defined by a boolean array where the **Cells** have been flagged by another **Filter**.  Typically, during reading/processing of the data, each **Cell** is subject to a "quality metric" (or threshold) that defines if the **Cell** is *good*.  This threshold can be used to define areas of each slice that are bad, either due to actual **Features** in the microstructure or external references inserted by the user/experimentalist.  If these "regions" of *bad* **Cells** are believed to be consistent through sections, then this **Filter** will preserve that by aligning those "regions" on top of one another on consecutive sections. 

This initial alignment is followed by a correction algorithm which adapts the shifts to a complementary criterion having one of the following forms:

1. Images of the mapped area taken from each slice by scanning electron microscope (EDAX).
2. User-defined shifts between the first and the last slice of the stack.

The initial alignment algorithm of this **Filter** is as follows:

1. Compare the value of the boolean array for each **Cell** in a section with the value of the array for the **Cell**  directly above it in the next section. 
2. Count the number of **Cell** pairs that do not have the same value and store that as the misalignment value for that position.
3. Repeat steps 1 and 2 for each position when shifting the second slice (relative to the first) from three (3) **Cells** to the left to three (3) **Cells** to the right, as well as from three (3) **Cells** up to three (3) **Cells** down. *Note that this creates a 7x7 grid*.
4. Determine the position in the 7x7 grid that has the lowest misalignment value. (It will be the position with the fewest different **Cell** pairs).
5. Repeat steps 1-4 with the center of each (new) 7x7 grid at the best position from the last 7x7 grid until the best position in the current/new 7x7 grid is the same as the last 7x7 grid. 
6. Repeat steps 1-5 for each pair of neighboring sections.

**Note that this is similar to a downhill simplex and can get caught in a local minimum!**

The correction alignment algorithm of this **Filter** attempts to improve the complementary fit as follows:

1. Start with the shifts obtained by the initial algorithm, i.e., define the current shifts as those obtained from the initial algorithm for each pair of consecutive cross sections.
2. For each pair of consecutive cross sections, consider all the shifts checked within the alignment algorithm (7x7 grids), different from the current shift, which improve the complementary fit. From these shifts, find the candidate shift as that with the lowest misalignment value. 
3. From the candidate shifts identified in step 2 (one candidate shift for each pair of consecutive sections), select the one with minimum difference between the misalignment value of the candidate shift and the misalignment value of the current shift. Change the current state by this candidate shift and recompute the complementary fit.
4. Repeat steps 2 - 3 until no improvement in the complementary fit is attained in step 2.

The user can choose to write the determined shift to an output file by enabling *Write Alignment Shifts File* and providing a file path.  


## Parameters ##

| Name | Type | Description |
|------|------| ----------- |
| Write Alignment Shift File | bool | Whether to write the shifts applied to each section to a file |
| Alignment File | File Path | The output file path where the user would like the shifts applied to the section to be written. Only needed if *Write Alignment Shifts File* is checked |
| Global Correction: SEM Images | bool | Whether to use SEM images for adaptive alignment. |
| Global Correction: Own Shifts | bool | Whether to set the shifts for adaptive alignment manually. |
| Total Shift In X-Direction (Microns) | float | Shift in X-direction between the first and the last slice of the stack in microns. Only needed if *Global Correction: Own Shifts* is checked. |
| Total Shift In Y-Direction (Microns) | float | Shift in Y-direction between the first and the last slice of the stack in microns. Only needed if *Global Correction: Own Shifts* is checked. |

## Required Geometry ##

Image 

## Required Objects ##

| Kind | Default Name | Type | Component Dimensions | Description |
|------|--------------|------|----------------------|-------------|
| **Image Data Array** | ImageData | uint8_t | any component size | Specifies image data containing the SEM images to be used for global correction of the shifts. |
| **Cell Attribute Array** | Mask | bool | (1) | Specifies if the **Cell** is to be counted as part of the *sample* or not. |

## Created Objects ##

None

## Example Pipelines ##



## License & Copyright ##

Please see the description file distributed with this **Plugin**

## DREAM.3D Mailing Lists ##

If you need more help with a **Filter**, please consider asking your question on the [DREAM.3D Users Google group!](https://groups.google.com/forum/?hl=en#!forum/dream3d-users)


