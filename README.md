# 4polar3D_SMOLM


# 🚧 Repository Status: Active Updates in Progress 🚧

> This repository is currently under **active development and code cleaning.**  
> Some scripts and workflows may still be updated or reorganized.  


## Overview  
**4polar3D_SMOLM** is a MATLAB-based toolkit for the analysis of **4polar3D microscopy datasets using polarized projection imaging.**  
It enables **preprocessing, polarized projection alignment, and trajectory coupling** to study **molecular orientation in complex biological samples.**

---

## Preliminary Setup Steps  
Before the main analysis begins, you must run the following scripts to prepare and correct the raw TIFF data:

1. `A_Make_polar_fbeads_stack.m`  
   Creates calibration stacks from polarized bead TIFF images.

2. `B_Rename_tiff_Hamamatsu.m`  
   Renames TIFF files acquired with Hamamatsu cameras to the correct stack format.

3. `C_Setregions.m`  
   Allows manual selection of the four polarized projection regions in the image.

4. `D_Correct_distortion.m`  
   Applies distortion correction using transformation matrices.

5. `E_Merge_tforms.m`  
   Merges multiple transformation files into a single transformation matrix for further processing.

6. `F_PRESSFIRST.m`  
   Sets up the directory structure and adds all necessary subfolders to the MATLAB search path.

---

## Processing Workflow
Once the data preparation steps are completed, follow this processing pipeline:

1. **Preprocess image stacks:**  
   - `FIRST_4polar3D_Preprocess.m`
   
2. **Analyze polarized projection pairs:**  
   - `SECOND_4polar3D_AnalyzePair_0_90_MultiCells.m` (0°/90°)  
   - `THIRD_4polar3D_AnalyzePair_45_135_MultiCells.m` (45°/135°)  

3. **Couple the projection pairs:**  
   - `FOURTH_4polar3D_CouplePairs_MultiCells.m`

---

## Key Features
- 📂 **Multi-stack processing:** Batch analysis across multiple cells or datasets.  
- 🔧 **Polarization-specific alignment:** Scripts for 0°/90° and 45°/135° projection pairs.  
- 🔁 **Distortion correction:** Uses transformation matrices to correct image distortions between projections.  
- 🔗 **Trajectory coupling:** Merges and analyzes vector trajectories across polarized datasets.

---

## Requirements
- MATLAB (recommended version: R2020a or newer)  
- Image Processing Toolbox  
- Parallel Computing Toolbox (recommended for preprocessing)

---

## Usage Notes
- ROIs must be carefully selected and dragged across paired polarized projections to ensure identical pixel coverage.
- Proper configuration of `STHERM_param.m` is required before running any analysis scripts.
- The coupling step matches trajectories between corresponding polarized projections.
- The pipeline is designed for multi-cell batch analysis.

---

## License  
This project is licensed under the **BSD 3-Clause "New" or "Revised" License.**
