# FTIR-Spectral-Analysis

Updated 13 September 2023
Author: Macy Mora-Antoinette

MATLAB R2023a  

<img width="844" alt="ftir" src="https://github.com/macy-mora-antoinette/FTIR-Spectral-Analysis/assets/112992304/f311315c-9bcd-4f10-93b5-8b8ea66df969">

This matlab code is intended to help users compute the compositional
properties from FTIR Absorption spectra in bones. 

Reference:
Hunt, H. B. et al. Altered Tissue Composition, Microarchitecture, and 
Mechanical Performance in Cancellous Bone From Men With Type 2 Diabetes 
Mellitus. Journal of Bone and Mineral Research 34, 1191–1206 (2019).
https://doi.org/10.1002/jbmr.3711

It is separated into three parts:

PART 1: A simple for loop to iterate through raw data files (csv)
in a single directory and take the average of all replicates. 
Plots are produced to visualize the average. 

PART 2: Processing the data using the average. Computes the area and peak
values for the different compositional properties. It will save plots for each
area and peak evaluated. Plots are to help the user evaluate data quality.  

PART 3: These functions used in PART 2 to help compute the different
values: spectra_peaks, spectra_area, and baseline.

Key things to keep in mind:
1. Works best if you run in sections PART 1 followed by PART 2
2. This code requires a specific naming format for your csv files. It
expects the sample replicates to have identical names except for #.csv at
the end to indicate the replicate.
Ex filename1.csv, filename2.csv, filename3.csv
2. It does not deal with corner cases, so you may need to manually change
wavelengths inputs for areas/peaks as needed
3. If you want to add additional area or peak computations, make sure you
change the "outputs" variable on line 134 and add the necessary columns to
the table on line 143
