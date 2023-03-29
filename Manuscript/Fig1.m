% wavelet coherence 
% This code is used to produce the figures included in 'wavelet coherence' section of Chapter 2 
% This code will not be submitted with the thesis 

% 
hlines = 14; 
inDir = '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Code (Github)/Icarus/inFiles/';
inDataFlame = '//Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Data/spectra/'; % Directory to spectra recorded with FLAME spectrometer at EOS
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

% Figure 1 
% a
% Intensity of two measured sky spectra at zenith 
inSpectra1 = '90_0_FLMS195681__0__10-13-49-165.txt'; % Spectrum 1
inSpectra2 = '90_0_FLMS195681__1__10-13-49-241.txt'; % Spectrum 2
s1 = fullfile(inDataFlame, inSpectra1); s1 = char(s1); s2 = fullfile(inDataFlame, inSpectra2); s2 = char(s2); 
if hlines >= 1
   [~, s1i] = textread(s1,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, s2i] = textread(s2,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [~, s1i] = textread(s1,'%f %f'); 
   [~, s2i] = textread(s2,'%f %f'); 
end  
figure
subplot(4,2,1) 
dp = plot(wHg, s1i, '-r'); hold on; dp.LineWidth = 1; % dp.Color = [0.8, 0, 0];
dp = plot(wHg, s2i, '-k'); hold on; dp.LineWidth = 1; % dp.Color = [0.8, 0, 0];
ylim([0 6e4]); ylabel('Intensity (counts)');
xlim([280 max(wHg)]); xlabel('\lambda (nm)'); set(gca,'XMinorTick','on','YMinorTick','on');
pos1 = [lAw1, 0, uAw1-lAw1, 6e4];
pos2 = [lAw2, 0, uAw2-lAw2, 6e4];
% rectangle('Position', pos1,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [0, 0, 0, 0.6] );
% rectangle('Position', pos2,'LineWidth', 0.5,'LineStyle', '--', 'EdgeColor', [0, 0, 0, 1] ); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% b
% Zoom in 
subplot(4,2,2)
dp = plot(wHg, s1i, '-r'); hold on; dp.LineWidth = 1; % dp.Color = [0.8, 0, 0];
dp = plot(wHg, s2i, '-k'); hold on; dp.LineWidth = 1; % dp.Color = [0.8, 0, 0];
pos1z = [lAw1, 1e4, uAw1-lAw1, 1e4];
pos2z = [lAw2, 1e4, uAw2-lAw2, 1e4];
% rectangle('Position', pos1z,'LineWidth', 1,'LineStyle', '-', 'EdgeColor', [0, 0, 0, 0.6] );
% rectangle('Position', pos2z,'LineWidth', 1,'LineStyle', '--', 'EdgeColor', [0, 0, 0, 0.4] ); 
ylim([1e4 2e4]); ylabel('Intensity (counts)'); set(gca, 'YScale', 'log')
xlim([lAw1 uAw1]); xlabel('\lambda (nm)'); set(gca,'XMinorTick','on'); 
Fig = gca; Fig.FontSize = 13; set(gcf,'color','w');
% c 
% Magnitude-squared wavelet coherence 
subplot(4,2,[3 6])
% wcoh = wcoherence(x,y) returns the magnitude-squared wavelet coherence, which is a measure of the correlation between signals x and y in the time-frequency plane. 
% The inputs x and y must be equal length, 1-D, real-valued signals. The coherence is computed using the analytic Morlet wavelet.
% [wcoh, wcs] = wcoherence(x,y) returns the wavelet cross-spectrum of x and y. 
% You can use the phase of the wavelet cross-spectrum values to identify the relative lag between the input signals.
% [wcoh, wcs,f] = wcoherence(x,y,fs) uses the positive sampling frequency, fs, to compute the scale-to-frequency conversion, f. 
% The sampling frequency fs is in Hz.
% [wcoh,wcs,f,coi] = wcoherence(___) returns the cone of influence, coi, for the wavelet coherence in cycles per sample. If you specify the sampling frequency, fs, the cone of influence is in Hz.
[wcoh, ~, fwc, coiwc] = wcoherence(s1i, s2i, Fs);
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([280 max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 16, 'Rotation', -90); 
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% Adds wavelength ranges of analysis window
Fh = max(fwc)-min(fwc); 
pos1wc = [lAw1, min(fwc), uAw1-lAw1, Fh];
pos2wc = [lAw2, min(fwc), uAw2-lAw2, Fh];
% rectangle('Position', pos1wc,'LineWidth', 1,'LineStyle', '-', 'EdgeColor', [0, 0, 0, 0.6] );
% rectangle('Position', pos2wc,'LineWidth', 1,'LineStyle', '--', 'EdgeColor', [0, 0, 0, 0.4] ); 

% % Adds window which will be extracted for later figures  
% maxF = 0.004;
% minF = 0.001;
% minW = 310;
% maxW = 330; 
% 
% % DEFINES ANALYSIS WINDOW
% % WAVELENGTH 
% % Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
% [~, indxEWl] = min(abs(wHg - minW)); % Lower 
% WlAw = wHg(indxEWl); 
% [~, indxEWu] = min(abs(wHg - maxW)); % Upper
% WuAw = wHg(indxEWu); 
% % FREQUENCY
% % Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
% [~, indxEFh] = min(abs(fwc - maxF)); % Highest 
% F1Aw = fwc(indxEFh); 
% [~, indxEFl] = min(abs(fwc - minF)); % Lowest 
% F2Aw = fwc(indxEFl); 
% 
% figure 
% [wcoh, ~, ~, ~] = wcoherence(s1i, s2i, Fs);
% wcoh(idxCoi) = NaN; 
% pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
% xlim([280 max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log'); % ylim([F2Aw F1Aw]); 
% cb = colorbar; 
% cb.Limits = [0 1];
% ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 16, 'Rotation', -90); 
% cb.Label.Position(1) = 4;
% Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
% 
% % Adds wavelength ranges of analysis window
% Fh = maxF-minF; 
% posE = [minW, minF, maxW-minW, Fh];
% rectangle('Position', posE,'LineWidth', 1,'LineStyle', '-', 'EdgeColor', [1, 1, 1] );






