# curvelet-denoising
Necessary functions and example script to perform curvelet denoising in MATLAB

# Inventory of scripts:
ExampleScript.m - script that integrates all scripts and items in this repository to perform an example of curvelet denoising of DAS data
CurveletDenoising.m - script that takes parameters for stochastic and coherent denoising to actually perform the denoising
MedianFilter.m - Script that performs median filtering as a preprocessing step

# Inventory of items:
EQ_raw.mat - MAT file of a DAS recording of an earthquake used in the example script
Noise_Matrix.mat - MAT file of a quiet DAS record section from the same array that recorded the earthquake in EQ_raw.mat
