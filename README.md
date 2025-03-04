# Histology-NMJ-and-muscle-Activity-Toolbox
A set of Matlab codes to:
  - Reconstruct Muscle 3D structures from muscle histology slides in nd2/tiff.
  - Overlaps features listed in excel file (this program was tested using NMJ locations).
  - Estimate volume and longest side of the muscle.
  - Show distribution patterns of features in excel using voxel grid and density based analysis.
  - Overlap high density EMG (hdEMG) activity recorded from muscle with histology reconstruction using signal attenuation modelling.
  - Correlate hdEMG activity to feature location using a electrode based signal reconstruction.
  - Analyzing neuromuscular junction (NMJ) distributions, generating EMG heatmaps, and visualizing electrode placements within muscle structures.

## Installation

To install and use these codes in MATLAB you can either clone the repository or download and access them.

# Main codes
## analysis_nd2.m and analysis_tiff.m
Analyses the histology slides in either tif or nd2 format and scales them properly, producing a 3D visualization of the slides stacked. It then plots and color codes the NMJ based on an excel sheet present in the same folder. 

### Requirements
- CY Y (2025). nd2read (https://github.com/joe-of-all-trades/nd2read), GitHub. Retrieved February 12, 2025.
- Aaron Wetzler (2025). Bresenham optimized for Matlab (https://www.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab), MATLAB Central File Exchange. Retrieved February 12, 2025.

## Avg_NMJ_contour.m
This code extracts and normalizes 3D point data from multiple muscle model figures, mapping them onto a common coordinate system based on a reference muscle model. It then defines a 3D voxel grid to calculate point densities and visualizes these densities using a heatmap overlay, along with plotting average point locations per voxel. Additionally, the script applies DBSCAN clustering to the normalized points, computes convex hulls for each dense cluster, and overlays these clusters onto the muscle model with transparent contour surfaces. Overall, it provides a detailed 3D visualization of muscle structure and NMJ distributions by integrating image data, point clustering, and spatial normalization.

### Requirements
- Oliver Woodford (2025). vol3d v2 (https://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2), MATLAB Central File Exchange. Retrieved February 12, 2025.
- extract3DPointsFromFig.m

## MeshOverlap.m
This code processes a 3D muscle reconstruction by first extracting NMJ (neuromuscular junction) coordinates and muscle structural data from a specified .fig file, and then loads a corresponding 2D heatmap that represents electrode activity. It constructs a mask to define an electrode region based on given scaling and positional parameters, overlays the scaled heatmap onto each muscle layer, and applies a centered subtraction to gradually decrease the intensity of already drawn areas. Finally, the code generates a 3D visualization by creating surfaces for each muscle slice with the adjusted electrode data and plots the NMJ coordinates on top.

### Requirements
- extract3DPointsFromFig.m

# Additional Features
## Muscle3DModel.m
Generates the convex hull of the 3D muscle reconstruction

### Requirements
- John Iversen (2025). freezeColors / unfreezeColors (https://github.com/jiversen/freezeColors/releases/tag/v2.5), GitHub. Retrieved February 17, 2025.

## ShowSide.m
Outlines the borders of the slides in the reconstruction to make them visible even from the side

## scalebar3D.m
Generates a 3D scalebar on a .fig file at coordinates [xc,yc,zc] with length w (number of units) and label.

### **CONE_INV_Dist_NMJ.m**
**Description:** Computes inverse distance mapping for NMJs in a conical muscle structure.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  
- Distance CSV file from `Distance_NMJ_CSV_Gen.m`  

### **CONE_VISUALISED.m**
**Description:** Visualizes the cone-shaped muscle model with NMJs and electrodes.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  

### **CYLINDER_INV_Dist_NMJ.m**
**Description:** Computes inverse distance mapping for NMJs in a cylindrical muscle structure.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  
- Distance CSV file from `Distance_NMJ_CSV_Gen.m`  

### **Distance_NMJ_CSV_Gen.m**
**Description:** Generates CSV files containing NMJ distance calculations.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  

### **HEATMAP_MUSCLE_OVERLAY.m**
**Description:** Overlays a generated EMG heatmap onto the muscle structure.  
**Requirements:**  
- Original muscle fig file with NMJs  
- Raw summed EMG heatmap from `RAW_HEATMAP_COMBINED.m`  

### **Heatmap_NMJ_Overlay_90_Shift.m**
**Description:** Creates an overlay of NMJs on an EMG heatmap with a 90-degree shift.  
**Requirements:**  
- Original muscle fig file with NMJs  
- Raw summed EMG heatmap from `RAW_HEATMAP_COMBINED.m`  

### **LOSSY_CONE_NMJ_Electrode_Point.m**
**Description:** Calculates electrode points and NMJ distances for a lossy cone model.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  

### **LOSSY_CONE.m**
**Description:** Implements a lossy cone-based model for NMJ distance computation.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  
- Distance CSV file from `LOSSY_CONE_NMJ_Electrode_Point.m`  

### **Oblate_Spheroid_INV_Dist_NMJ.m**
**Description:** Computes inverse distance mapping for NMJs in an oblate spheroid muscle model.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  
- Distance CSV file from `Distance_NMJ_CSV_Gen.m`  

### **OBLATE_SPHEROID_VISUALISED.m**
**Description:** Visualizes the oblate spheroid muscle model with NMJs and electrodes.  
**Requirements:**  
- Muscle fig file with NMJs and electrode positions from `Heatmap_NMJ_Overlay_90_Shift.m`  
