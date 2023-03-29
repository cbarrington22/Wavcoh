# Wavcoh

This repository contains the MATLAB code used to obtain the results presented in “Barrington et al. Simplifying the Retrieval of Volcanic SO2 from UV spectra using Wavelet Coherence”.
The results presented in this manuscript were produced using the various “.m” files contained in the ”Manuscript" folder which have been uploaded here in the interest of research integrity. This version of the code is not thoroughly commented nor is it free of redundancies. Figures produced using this code have been modified and formatted independently although the data remains unchanged.
Users wishing to apply similar analysis to their own spectra, are instead directed to the USER-FRIENDLY SCRIPT ("UserFriendly.m") which provides the necessary code to apply the method to UV spectra recorded by scanning UV instruments. This script calculates the magnitude-squared wavelet coherence (MSWC) between spectra to provide a measure of correlation in the wavelength-spatial frequency plane. The minimum MSWC and the differential MSWC are determined within an extraction window between spatial frequencies of approximately 1 to 4 · 10-3 cycles/lambda and wavelengths of approximately 310.0 and 326.8 nm, where the characteristics of the UV spectra are highly uniform but where the changes in the MSWC due to the narrow band absorption of SO2 are visible. The minimum MSWC and the differential MSWC are plotted against the scan elevation angle, also indicating the reference spectrum (please see Figures 7 and 9 in the manuscript for examples). The MSWC may also be plotted (similar to Figures 8, S2, S3 and S4 in the manuscript). The phase angle of the wavelet cross spectrum (WCS) is also determined (as in Figures S7 and S11) and the mean phase angle within the extraction window is plotted against the scan elevation angle. By default, this code will run with spectra recorded by the NOVAC (Network of Volcanic and Atmospheric Change) but may be adapted to fit the user’s needs. The code in this script differs slightly from that used to produce the figures included in the manuscript. Hardcoding and redundant lines have been removed (as far as possible), the analysis procedure is the same however. 

Code is written in MATLAB using version R2021a (MathsWorks, 2019) and utilises several in-built functions including the ‘wcoherence’ function to determine the MSWC and WCS. This function requires the Wavelet Toolbox. For further details see: https://www.mathworks.com/help/wavelet/ref/wcoherence.html.
