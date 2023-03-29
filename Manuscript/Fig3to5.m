% wavelet coherence 
% This code is used to produce the figures included in 'wavelet coherence' section of Chapter 2 
% This code will not be submitted with the thesis 

% 
hlines = 14; 
inDir = '/Users/charlotteb/Documents/Frontiers SIRS manuscript submission/Code (Github)/Icarus/inFiles/';
inDataFlame = '/Users/charlotteb/Documents/Frontiers SIRS manuscript submission/Data/spectra/'; % Directory to spectra recorded with FLAME spectrometer at EOS
% Primary and secondary wavelength ranges 
lAw1 = 314.8;
uAw1 = 326.8; 
lAw2 = 326.8; 
uAw2 = 353.5; 

% Loads Hg spectrum to obtain wavelet information consistent with that used in linear model  
fnHg = 'FLMS195681__0__12-14-16-726.txt'; % Hg spectrum 
dHg = fullfile(inDir, fnHg); inHg = char(dHg); % Hg spectrum 
if hlines >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  

% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequency is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

% Figure 3
% a
% Intensity of two measured sky spectra at zenith 
inSpectraS = '90_0_FLMS195681__0__10-13-49-165.txt'; % Sky 
inSpectraSs = '90_0_FLMS195681__1__10-13-49-241.txt'; % Sky 2
inSpectra1 = '90_51_FLMS195681__0__10-17-52-133.txt'; % Spectrum 1
inSpectra2 = '90_516_FLMS195681__0__10-20-00-029.txt'; % Spectrum 2
inSpectra3 = '90_1063_FLMS195681__0__10-21-30-626.txt'; % Spectrum 3
inSpectra4 = '90_2080_FLMS195681__0__10-16-37-635.txt'; % Spectrum 4
inSpectraA = '90_113_FLMS195681__0__10-20-51-127.txt'; % Additional
%
inSpectraA2 = '90_107_FLMS195681__0__10-19-23-430.txt'; % 
inSpectraA3 = '90_1019_FLMS195681__0__10-18-30-331.txt'; % 
inSpectraA4 = '90_1526_FLMS195681__0__10-17-16-234.txt'; % 
inSpectraA5 = '90_1937_FLMS195681__0__10-22-07-025.txt'; % 
% 
sS = fullfile(inDataFlame, inSpectraS); sS = char(sS); sSs = fullfile(inDataFlame, inSpectraSs); sSs = char(sSs);
s1 = fullfile(inDataFlame, inSpectra1); s1 = char(s1); s2 = fullfile(inDataFlame, inSpectra2); s2 = char(s2); 
s3 = fullfile(inDataFlame, inSpectra3); s3 = char(s3); s4 = fullfile(inDataFlame, inSpectra4); s4 = char(s4); 
sA = fullfile(inDataFlame, inSpectraA); sA = char(sA); 
% 
sA2 = fullfile(inDataFlame, inSpectraA2); sA2 = char(sA2); 
sA3 = fullfile(inDataFlame, inSpectraA3); sA3 = char(sA3); 
sA4 = fullfile(inDataFlame, inSpectraA4); sA4 = char(sA4); 
sA5 = fullfile(inDataFlame, inSpectraA5); sA5 = char(sA5); 
if hlines >= 1
   [~, sSi] = textread(sS,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sSsi] = textread(sSs,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, s1i] = textread(s1,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, s2i] = textread(s2,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, s3i] = textread(s3,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, s4i] = textread(s4,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sAi] = textread(sA,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   % 
   [~, sA2i] = textread(sA2,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sA3i] = textread(sA3,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sA4i] = textread(sA4,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sA5i] = textread(sA5,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [~, sSi] = textread(sS,'%f %f'); 
   [~, sSsi] = textread(sSs,'%f %f'); 
   [~, s1i] = textread(s1,'%f %f'); 
   [~, s2i] = textread(s2,'%f %f'); 
   [~, s3i] = textread(s3,'%f %f'); 
   [~, s4i] = textread(s4,'%f %f'); 
   [~, sAi] = textread(sA,'%f %f'); 
   % 
   [~, sA2] = textread(sA2,'%f %f'); 
   [~, sA3] = textread(sA3,'%f %f'); 
   [~, sA4] = textread(sA4,'%f %f'); 
   [~, sA5] = textread(sA5,'%f %f'); 
end  
% Magnitude-squared wavelet coherence 
figure
% wcoh = wcoherence(x,y) returns the magnitude-squared wavelet coherence, which is a measure of the correlation between signals x and y in the time-frequency plane. 
% The inputs x and y must be equal length, 1-D, real-valued signals. The coherence is computed using the analytic Morlet wavelet.
% [wcoh, wcs] = wcoherence(x,y) returns the wavelet cross-spectrum of x and y. 
% You can use the phase of the wavelet cross-spectrum values to identify the relative lag between the input signals.
% [wcoh, wcs,f] = wcoherence(x,y,fs) uses the positive sampling frequency, fs, to compute the scale-to-frequency conversion, f. 
% The sampling frequency fs is in Hz.
% [wcoh,wcs,f,coi] = wcoherence(___) returns the cone of influence, coi, for the wavelet coherence in cycles per sample. If you specify the sampling frequency, fs, the cone of influence is in Hz.
subplot(2,2,1) 
[wcoh, ~, fwc, coiwc] = wcoherence(sSi, sSsi, Fs);
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); hold on;
xlim([280 max(wHg)]); set(gca,'XMinorTick','on','YMinorTick','on'); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 12, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
title('No gas cell')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
% Adds window which will be extracted for figure 4 
maxF = 0.004;
minF = 0.001;
minW = 310;
maxW = 326.8; 
% DEFINES ANALYSIS WINDOW
% WAVELENGTH 
% Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
[~, indxEWl] = min(abs(wHg - minW)); % Lower 
WlAw = wHg(indxEWl); 
[~, indxEWu] = min(abs(wHg - maxW)); % Upper
WuAw = wHg(indxEWu); 
% FREQUENCY
% Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
[~, indxEFh] = min(abs(fwc - maxF)); % Highest 
F1Aw = fwc(indxEFh); 
[~, indxEFl] = min(abs(fwc - minF)); % Lowest 
F2Aw = fwc(indxEFl); 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
fwcEx = fwc(indxEFh:indxEFl, :); 
wHgEx = wHg(indxEWl:indxEWu, :); 
minMSWC = zeros(4, 1); 
minMSWC(1,1) = min(wcohEx,[],'all'); % Returns the smallest element
% Adds wavelength ranges of analysis window
Fh = maxF-minF; 
posE = [minW, minF, maxW-minW, Fh];
rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );

subplot(2,2,2) 
[wcoh, ~, fwc, ~] = wcoherence(sSi, s2i, Fs);
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); hold on;
xlim([280 max(wHg)]); set(gca,'XMinorTick','on','YMinorTick','on'); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 12, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
title('516 ppm·m')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
minMSWC(2,1) = min(wcohEx,[],'all'); % Returns the smallest element
% Adds wavelength ranges of analysis window
Fh = maxF-minF; 
posE = [minW, minF, maxW-minW, Fh];
rectangle('Position', posE,'LineWidth', 0.1,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );

subplot(2,2,3) 
[wcoh, ~, fwc, ~] = wcoherence(sSi, s3i, Fs);
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); hold on;
xlim([280 max(wHg)]); set(gca,'XMinorTick','on','YMinorTick','on'); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 12, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
title('1063 ppm·m')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
minMSWC(3,1) = min(wcohEx,[],'all'); % Returns the smallest element
% Adds wavelength ranges of analysis window
Fh = maxF-minF; 
posE = [minW, minF, maxW-minW, Fh];
rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );

subplot(2,2,4) 
[wcoh, ~, fwc, ~] = wcoherence(sSi, s4i, Fs);
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); hold on;
xlim([280 max(wHg)]); set(gca,'XMinorTick','on','YMinorTick','on'); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 12, 'Rotation', -90); 
cb.Label.Position(1) = 3.5;
title('2080 ppm·m')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
minMSWC(4,1) = min(wcohEx,[],'all'); % Returns the smallest element
% Adds wavelength ranges of analysis window
Fh = maxF-minF; 
posE = [minW, minF, maxW-minW, Fh];
rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );

% Figure 4 
% Zoom in - extraction window
figure
% wcoh = wcoherence(x,y) returns the magnitude-squared wavelet coherence, which is a measure of the correlation between signals x and y in the time-frequency plane. 
% The inputs x and y must be equal length, 1-D, real-valued signals. The coherence is computed using the analytic Morlet wavelet.
% [wcoh, wcs] = wcoherence(x,y) returns the wavelet cross-spectrum of x and y. 
% You can use the phase of the wavelet cross-spectrum values to identify the relative lag between the input signals.
% [wcoh, wcs,f] = wcoherence(x,y,fs) uses the positive sampling frequency, fs, to compute the scale-to-frequency conversion, f. 
% The sampling frequency fs is in Hz.
% [wcoh,wcs,f,coi] = wcoherence(___) returns the cone of influence, coi, for the wavelet coherence in cycles per sample. If you specify the sampling frequency, fs, the cone of influence is in Hz.
subplot(4,3,1) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sSsi, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
MEAN = mean(wcohEx, 'all');
SD = std(wcohEx(:));
MAX = max(wcohEx(:));
SUM = sum(wcohEx(:));

pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('No gas cell')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx = zeros(10, 1); 
totMSWCEx(1,1) = sum(wcohEx,'all'); 
MINtotMSWCEx = zeros(10, 1); 
SDtotMSWCEx = zeros(10, 1); 
MINtotMSWCEx(1,1) = min(wcohEx(:)); 
SDtotMSWCEx(1,1) = std(wcohEx(:)); 

subplot(4,3,2) 
[wcoh, ~, ~, ~] = wcoherence(sSi, s1i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('51 ppm·m')
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
totMSWCEx(2,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(2,1) = min(wcohEx(:)); 
SDtotMSWCEx(2,1) = std(wcohEx(:)); 

subplot(4,3,3) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sA2i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('107 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(3,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(3,1) = min(wcohEx(:)); 
SDtotMSWCEx(3,1) = std(wcohEx(:)); 

subplot(4,3,4) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sAi, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('113 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(4,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(4,1) = min(wcohEx(:)); 
SDtotMSWCEx(4,1) = std(wcohEx(:)); 

subplot(4,3,5) 
[wcoh, ~, ~, ~] = wcoherence(sSi, s2i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('516 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(5,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(5,1) = min(wcohEx(:)); 
SDtotMSWCEx(5,1) = std(wcohEx(:)); 

subplot(4,3,6) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sA3i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('1019 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(6,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(6,1) = min(wcohEx(:)); 
SDtotMSWCEx(6,1) = std(wcohEx(:)); 

subplot(4,3,7) 
[wcoh, ~, ~, ~] = wcoherence(sSi, s3i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('1063 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(7,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(7,1) = min(wcohEx(:)); 
SDtotMSWCEx(7,1) = std(wcohEx(:)); 

subplot(4,3,8) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sA4i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]);set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('1526 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(8,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(8,1) = min(wcohEx(:)); 
SDtotMSWCEx(8,1) = std(wcohEx(:)); 

subplot(4,3,9) 
[wcoh, ~, ~, ~] = wcoherence(sSi, sA5i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('1937 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(9,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(9,1) = min(wcohEx(:)); 
SDtotMSWCEx(9,1) = std(wcohEx(:)); 

subplot(4,3,10) 
[wcoh, ~, ~, ~] = wcoherence(sSi, s4i, Fs);
wcoh(idxCoi) = NaN; 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]); set(gca,'YScale','log'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'MSWC', 'FontSize', 11, 'Rotation', -90); 
cb.Label.Position(1) = 4;
title('2080 ppm·m')
Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
totMSWCEx(10,1) = sum(wcohEx,'all'); 
MINtotMSWCEx(10,1) = min(wcohEx(:)); 
SDtotMSWCEx(10,1) = std(wcohEx(:)); 

totMSWCEx = round(totMSWCEx,2,"decimals");
disp(totMSWCEx)

% Figure 5 
% Differencial coherence 
blMSWC = totMSWCEx(1,1); 
d = totMSWCEx-blMSWC;

% Gas cell concentration 
gacConAll = [0; 51; 107; 113; 516; 1019; 1063; 1526; 1937; 2080]; 
gasCon1 = [113; 1063; 2080]; 
gasCon2 = [0; 51; 107; 516; 1019; 1526; 1937]; 
d1 = [d(4); d(7); d(10)];
d2 = [d(1); d(2); d(3); d(5); d(6); d(8); d(9)];
% Find 10% of gasCon 
eR = (gacConAll/100)*10; 
figure 
subplot(2,1,2)
p = plot(gacConAll(1:end,:), abs(d(1:end,:)), '.k'); hold on;
p.MarkerSize = 8; 
h = lsline; % Superimposes a least-squares line on each scatter plot in the current axes.
h.LineWidth = 1;
% Pause and delete p 
p = errorbar(gacConAll(2:end,:),abs(d(2:end,:)), eR(2:end,:), eR(2:end,:), 'horizontal', '.k'); hold on;
p.CapSize = 0; 
p.Marker = 'none'; 
p.LineStyle = 'none'; 
p = plot(gasCon2(1:end,:), abs(d2(1:end,:)), '.k'); hold on;
p.MarkerSize = 8; 
p = plot(gasCon1(1:end,:), abs(d1(1:end,:)), '.k'); hold on;
p.MarkerSize = 8; 
xlabel('Gas cell concentration (ppm·m)'); ylabel('Differential MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');


% Figure 5 
% Differencial coherence 
blMSWC = totMSWCEx(1,1); 
d = 1-MINtotMSWCEx;
% dErr = SDtotMSWCEx;
subplot(2,1,1)
% Gas cell concentration 
gacConAll = [0; 51; 107; 113; 516; 1019; 1063; 1526; 1937; 2080]; 
gasCon1 = [113; 1063; 2080]; 
gasCon2 = [0; 51; 107; 516; 1019; 1526; 1937]; 
d1 = [d(4); d(7); d(10)];
d2 = [d(1); d(2); d(3); d(5); d(6); d(8); d(9)];
% Find 10% of gasCon 
eR = (gacConAll/100)*10; 
% p = errorbar(gacConAll(1:end,:), abs(d(1:end,:)), abs(dErr(1:end,:)), abs(dErr(1:end,:))); hold on 
% eR(2:end,:), eR(2:end,:), 'k'); hold on;
p = plot(gacConAll(1:end,:), abs(d(1:end,:)), '.k'); hold on;
% h = lsline; % Superimposes a least-squares line on each scatter plot in the current axes.
% h.LineWidth = 1;
% Pause and delete p 
p = errorbar(gacConAll(2:end,:),abs(d(2:end,:)), eR(2:end,:), eR(2:end,:), 'horizontal', '.k'); hold on;
p.CapSize = 0; 
p.Marker = 'none'; 
p.LineStyle = 'none'; 
p = plot(gasCon2(1:end,:), abs(d2(1:end,:)), '.k'); hold on;
p.MarkerSize = 8; 
p = plot(gasCon1(1:end,:), abs(d1(1:end,:)), '.k'); hold on;
p.MarkerSize = 8; 
% xlabel('Gas cell concentration (ppm·m)'); ylabel('1/min. MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');



