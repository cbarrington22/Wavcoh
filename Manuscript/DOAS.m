% Loads DOAS data 
addpath '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'
inDirD = '/Users/charlotteb/Documents/Chapter 2/NOVAC spectra/'; 
inFileD_1417 = 'mDOASResults_2022-09-27T11-37-44_D2J2375_140319_1417_0.csv'; 
inFileD_1941 = 'mDOASResults_2022-09-27T11-40-34_D2J2375_140319_1941_0.csv'; 
hlinesD = 1; 
% 1417 
inD = fullfile(inDirD, inFileD_1417); 
% Loads all data in results file
[SpectrumID_1417, FitCoeff_6101_1417,	FitCoeffError_6101_1417,	Shift_6101_1417,	Squeeze_6101_1417, FitCoeff_Ring_1417, FitCoeffError_Ring_1417,	Shift_Ring_1417,...
    Squeeze_Ring_1417, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1417, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417, nSO2D_1417, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417,...
    WavelengthRangeLow_1417,	WavelengthRangeHigh_1417, FitChi2_1417, TimeStamp_1417, Latitude_1417, Longitude_1417, Altitude_1417,	Speed_1417, Course_1417, GPSWarnCode_1417,	GPSQuality_1417,	ElevationAngle_1417,	AzimuthAngle_1417,...
    ExposureTime_1417, Exposures_1417, MaxIntensity_1417, MaxIntensityFitRange_1417, FileName_1417, SpectrometerType_1417, SpectrometerSerialNumber_1417,...
    SpectrometerChannel_1417, Remark_1417] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
% 1941 
inD = fullfile(inDirD, inFileD_1941); 
[SpectrumID_1941, FitCoeff_6101_1941,	FitCoeffError_6101_1941,	Shift_6101_1941,	Squeeze_6101_1941, FitCoeff_Ring_1941, FitCoeffError_Ring_1941,	Shift_Ring_1941,...
    Squeeze_Ring_1941, FitCoeff_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,	FitCoeffError_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, Shift_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941,...
    Squeeze_FLMS195681_O3_Voigt_223K_HP500_PPMM_1941, FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941,...
    FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941, nSO2D_1941, Shift_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,	Squeeze_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941,...
    WavelengthRangeLow_1941,	WavelengthRangeHigh_1941, FitChi2_1941, TimeStamp_1941, Latitude_1941, Longitude_1941, Altitude_1941,	Speed_1941, Course_1941, GPSWarnCode_1941,	GPSQuality_1941,	ElevationAngle_1941,	AzimuthAngle_1941,...
    ExposureTime_1941, Exposures_1941, MaxIntensity_1941, MaxIntensityFitRange_1941, FileName_1941, SpectrometerType_1941, SpectrometerSerialNumber_1941,...
    SpectrometerChannel_1941, Remark_1941] = textread(inD,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%f%s', 'headerlines', hlinesD, 'delimiter', ',');  
 
SO2D_1417 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1417;  
SO2ErrD_1417 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1417; 

SO2D_1941 = FitCoeff_FLMS195681_SO2_Bogumil_293K_HP500_MOLECCM2_1941;  
SO2ErrD_1941 = FitCoeffError_FLMS195681_SO2_Bogumil_293K_HP500_PPMM_1941; 

yyaxis right 
e = errorbar(ElevationAngle_1417(cSpec(3:end)), nSO2D_1417(cSpec(3:end)), SO2ErrD_1417(cSpec(3:end)), SO2ErrD_1417(cSpec(3:end)), '-ob'); hold on; % Plots minimum wavelet coherend against scan angle 
e.MarkerFaceColor = 'b'; 
e.MarkerEdgeColor = 'k'; 
ylabel('SCD SO_2 (molec/cm^2)'); 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend('1-min MSWC', 'SCD SO_2');

yyaxis right 
e = errorbar(ElevationAngle_1941(cSpec(3:end)), nSO2D_1941(cSpec(3:end)), SO2ErrD_1941(cSpec(3:end)), SO2ErrD_1941(cSpec(3:end)), '-ob'); hold on; % Plots minimum wavelet coherend against scan angle 
ylabel('SCD SO_2 (molec/cm^2)'); 
e.MarkerFaceColor = 'b'; 
e.MarkerEdgeColor = 'k'; 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend('1-min MSWC', 'SCD SO_2');
