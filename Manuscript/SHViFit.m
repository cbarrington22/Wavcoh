% wavelet coherence 
% This code is used to produce the figures included in 'wavelet coherence' section of Chapter 2 
% This code will not be submitted with the thesis 

% SHV scan
% Take one full instrument scan 
% Select the scan with XX degree elevation angle as the reference (informed from iFit analysis)
% Wavelet coherence of each subsequent scan (to 60 degrees) 
% We don't have Hg spectrum for SHV so we use channel 
% Extract analysis window accord to previously used window 
% What is the wavelet coherence threshold identified by the zenoth-60 degree elevation angle test 
% That's the baseline - we only consider wavelet coherence which crosses that threshold total wc 5.4215e+03 (Fig. 6 ) - assume same for this spectrometer 
% Plot a pannel of subplots showing SO2 absorption signature 
% Extract minimum total wavelet coherence of the extraction window 
% Save total wavelet coherence and elevtaion angle of each scan 

inDir = '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Code (Github)/Icarus/inFiles/';
toAnalyse = '/Users/charlotteb/Desktop/Thesis/Content/Chapter 2/Data/SHV/USB2+F03187_190516_1428_0/'; % Directory to full full instrument scan

% Copy of Appendix A 
fnSolar = 'SAO2010_Chance_Kurucz.txt'; % Solar reference
dSolar = fullfile(inDir, fnSolar); inSolar = char(dSolar); [wSolar, aSolar] = textread(inSolar,'%f %f'); % Solar reference

hlines = 14; 

% Loads Hg spectrum to obtain wavelet information consistent with that used in linear model  
fnHg = 'FLMS195681__0__12-14-16-726.txt'; % Hg spectrum 
dHg = fullfile(inDir, fnHg); inHg = char(dHg); % Hg spectrum 
if hlines >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end

% Copy of Appendix A 
fnSolar = 'SAO2010_Chance_Kurucz.txt'; % Solar reference
dSolar = fullfile(inDir, fnSolar); inSolar = char(dSolar); [wSolar, aSolar] = textread(inSolar,'%f %f'); % Solar reference

% % RESAMPLES DATA AT WAVELENGTHS DEFINED BY CALIBRATED Hg spectrum
% % Crops wHg to limit wavelengths to those 280 nm or above
% [ind, ~] = find(wHg >= 300);
% wHg = wHg(ind, 1); 
% iHg = iHg(ind, 1); 
 
% Uses interp(x, v, xq) to query trace gas cross sections and solar reference at wavelengths defined by calibrated Hg spectrum   
% 'interp' returns interpolated values of a 1-D function at specific query points (xq) using linear interpolation 
% Vector x contains the sample points - original wavelength information 
% Vector v contains the corresponding values, v(x) - original absorption (OR intensity) information 
solar = interp1(wSolar, aSolar, wHg,'linear'); 

% CONVOLUTION
% Convolutes resampled trace gas cross sections and solar reference with data kernel to match ILS of spectrometer 
figure('Renderer', 'painters', 'Position', [900 900 900 600]) % Plots Hg spectrum 
dp = plot(wHg, iHg); hold on; dp.LineWidth = 1;
xlim([min(wHg) max(wHg)]); xtickn = (280:10:420); xticks(xtickn);
set(gca,'XMinorTick','on'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)');
legend({'Recorded Hg spectrum','Wavelength range of analysis window'},'FontSize',14);  
title('Selecting Hg-peak to use for data kernel k(x)'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% Defines data kernel k(x)
pL = sprintf('Select a single peak of the Hg spectrum to be used for the data kernel.\n Select a peak which is not saturated but in close proximity to the analysis window.\n Input the x-axis value corresponding to the start of the selected Hg-peak (l):\n');
l = input(pL);
pH = sprintf('Input the x-axis value corresponding to the end of the selected Hg-peak (h):\n');
h = input(pH);
% dp = xline(lAw,'-k'); dp.LineWidth = 1;
dp = xline(l,'-r'); dp.LineWidth = 1;
% dp = xline(uAw, '-k'); dp.LineWidth = 1;
dp = xline(h, '-r'); dp.LineWidth = 1;

% Extracts data s(x) from k(x)
[rI, ~] = find(wHg > l&wHg < h); % Returns index of wavelengths between defined limits of kx (iL and iH)
sxW = wHg(rI, :); % Extracts wavelength data for kernel 
sxI = iHg(rI, :); % Extracts intensity data for kernel  
% Plots data kernel s(x) with Smin
figure('Renderer','painters','Position',[900 900 900 600]) 
subplot(1,3,1) % Plots data kernel s(x) with Smin as red cross
dp = plot(sxW, sxI); hold on; dp.LineWidth=1;
[smin, sminInd] = min(sxI);
plot(sxW(sminInd, :), smin, 'xr');
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
subplot(1,3,2) % Plots data kernel s(x) - Smin 
sxISmin = sxI - smin;
dp = plot(sxW, sxISmin); hold on; dp.LineWidth=1;
xlim([min(sxW) max(sxW)]);
xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); legend({'s(x) - S_m_i_n from l to h'},'FontSize', 14, 'Location', 'NorthOutside');    
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
subplot(1,3,3) % Plots normalised data kernel k(x) =  s(x) - Smin / ‚à´l^h s(x) - Smin
% 'trapz' computes the approximate integral of Y via the trapezoidal method with unit spacing. If Y is a vector, then trapz(Y) is the approximate integral of Y
kx = sxISmin/(trapz(sxISmin)); 
dp = plot(sxW,kx); hold on; dp.LineWidth = 1;
xlim([min(sxW) max(sxW)]); xlabel('\lambda (nm)'); ylabel('Normalised intensity');
legend({'s(x) - S_m_i_n / ‚à´_l^h s(x) - S_m_i_n'},'FontSize', 14,'Location', 'NorthOutside');   
t = sgtitle('Defining data kernel k(x) for convolution'); t.FontWeight = 'bold'; 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
 
% Convolutes resampled trace gas cross sections
% w = conv(u, v, shape) returns a subsection of the convolution, as specified by shape. For example, conv(u, v, 'same') returns only the central part of the convolution, the same size as u
cSolar = conv(solar, kx, 'same');
 
% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequency is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

% 
subFolders=dir(toAnalyse); % Returns data on subfolders to be analysed
subFoldersData=subFolders(~ismember({subFolders(:).name},{'.','..','.DS_Store'})); % Removes '.' and '..' and '.DS_Store' contained in subFolder varaible
spectraList=extractfield(subFoldersData,'name');spectraList=spectraList'; % Extracts file names to list all spectra to be analysed
dark=[];
sky=[];
scan=[];
angleOut=[];

% Loads wavelength information 
a = 1; 
spectraName=spectraList(a); 
spectraDirectory=fullfile(toAnalyse,spectraName); % Directory to access individual spectra within directory
inSpectraI=char(spectraDirectory);
[data]=textread(inSpectraI,'%s');
name=data(2145,:); % For some reason this is different % 2150
name=string(name);
angle=data(2091,:); 
spectra=data(4:2049,:); 
spectra=str2double(spectra);

% PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
colsolar = [0 0 0];
figure 
yyaxis left 
dp = plot(wHg, cSolar, '-'); hold on; dp.LineWidth = 1; dp.Color = colsolar;
xlim([310 325]);
ylim([0 max(cSolar)]);
ylabel('Intensity (counts)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
% dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
% dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
yyaxis right 
plot(wHg(1:end-4)+7, spectra(1:end-2)); hold on;
ylabel('Intensity (counts)');
legend('Solar reference', 'NOVAC scan')
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
dp = xline(310,'-r'); dp.LineWidth = 1;
dp = xline(325,'-r'); dp.LineWidth = 1;

for a=1:length(spectraList) 
    spectraName=spectraList(a); 
    spectraDirectory=fullfile(toAnalyse,spectraName); % Directory to access individual spectra within directory
    inSpectraI=char(spectraDirectory);
    [data]=textread(inSpectraI,'%s');
    name=data(2145,:); 
    name=string(name);
    angle=data(2091,:); 
    spectra=data(4:2047,:); % 2051-4
    spectra=str2double(spectra);
    switch name 
        case '"sky"'
        sky=[sky,spectra];
        case '"dark"'
        dark=[dark,spectra];
        case  '"scan"'
        scan=[scan,spectra];
        angleOut=[angleOut,angle];
        otherwise 
        warning('Unexpected plot type. No plot created.');
    end 
end 

% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequency is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg(1:end-4)); % Number of spectrometer channels  
minW = min(wHg(1:end-4)); % Starting wavelength 
maxW = max(wHg(1:end-4)); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

% convert angle to number 
scanAngle = str2double(angleOut);

% Reference angle 
% targrefAng = 60; 
% % [~,idx] = min(abs(scanAngle - targrefAng));
% % refAng = scanAngle(idx);
% % refAng = scanAngle(9);
idx = 9; 
refSpec = scan(:, idx); % Reference scan -60 
refAng = scanAngle(idx);

% [~,cSpec] = find((scanAngle>refAng) & scanAngle<(60)); Limitted upper scan angle  
[~,cSpec] = find(scanAngle>refAng);

[~, ~, fwc, ~] = wcoherence(refSpec, scan(:, idx+1), Fs);

run = 1;
maxF = 0.005; 
minF = 0.001;
minW = 310; % 300;
maxW = 325; % 326.8; 

wHg = wHg(1:end-4)+7; 

% Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
[~, indxEWl] = min(abs(wHg - minW)); % Lower 
WlAw = (indxEWl); 
[~, indxEWu] = min(abs(wHg - maxW)); % Upper
WuAw = wHg(indxEWu); 
% FREQUENCY
% Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
[~, indxEFh] = min(abs(fwc - maxF)); % Highest 
F1Aw = fwc(indxEFh); 
[~, indxEFl] = min(abs(fwc - minF)); % Lowest 
F2Aw = fwc(indxEFl); 

[wcoh, ~, fwc, coiwc] = wcoherence(refSpec, scan(:, idx+1), Fs);
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI

maxWC = 0.8; % 0.9871; % from figure 6 
out = []; 

figure 
for i = cSpec % 1:cSpec 
    y = scan(:, i);
    [wcoh, ~, ~, ~] = wcoherence(refSpec, y, Fs);
    wcoh(idxCoi) = NaN; 
    % subplot(7, 5, run)
    subplot(7, 6, run)
    pcolor(wHg, fwc, abs(wcoh)); shading interp; colormap(hot(50)); set(gca, 'YScale', 'log')
    xlim([300 326.8]); ylim([minF maxF]); 
    box on
    Ang = string(scanAngle(:, i));
    title(sprintf('%s°', Ang));
    wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
    minn = min(wcohEx, [], 'all');
    out = [out, minn]; 
    xticks([300 310 320])
    xticklabels({' ','310','320'}); set(gca,'XMinorTick','on')
    Fig = gca; Fig.FontSize = 7; set(gcf,'color','w');
    run = run + 1;  
end 

cb = colorbar; 
cb.Limits = [0 maxWC];

% % Figure 5 
% % Differencial coherence 
SHV = out;
% % masaya2 = out;
% 
figure % All in window 
x = scanAngle(:,cSpec); 
yyaxis left 
plot(x, 1-SHV, '-ok'); hold on; % Plots minimum wavelet coherend against scan angle 
% plot(x, 1-masaya2, '-ob'); hold on; % Plots minimum wavelet coherend against scan angle (plot ifit) 
xlabel('Scan angle (°)'); ylabel('1 - minimum MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

% Loads NOVAC and iFit data 
ifit = 'iFit_USB2+F03187_190516_1428_0_mod.txt'; 
inDir = '/Users/charlotteb/Desktop/Thesis/Content/Chapter 2/Data/SHV'; 

ifit = fullfile(inDir, ifit); ifit = char(ifit); % Hg spectrum 
hlines = 1; 
[Number, ifitSO2, ifitSO2_err] = textread(ifit, '%f%f%f', 'headerlines', hlines); % Skips header lines if exist 
yyaxis right
errorbar(scanAngle(:,cSpec), ifitSO2(cSpec,:), ifitSO2_err(cSpec,:), ifitSO2_err(cSpec,:), '-or'); hold on; % Plots minimum wavelet coherend against scan angle 
ylabel('SCD SO_2 (molec/cm^2)'); 
% Loads NOVAC and iFit data 
novac = 'NOVAC_USB2+F03187_190516_1428_0_mod.txt';
novac = fullfile(inDir, novac); novac = char(novac); % Hg spectrum 
[sAngle, SO2c, colErr] = textread(novac, '%f%f%f', 'headerlines', hlines); % Skips header lines if exist 
% errorbar(sAngle(cSpec+1,:), SO2c(cSpec+1,:), colErr(cSpec+1,:), colErr(cSpec+1,:), '-ob'); hold on; % Plots minimum wavelet coherend against scan angle 
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';




