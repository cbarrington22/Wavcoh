% wavelet coherence 
% This code is used to produce the figures included in 'wavelet coherence' section of Chapter 2 
% This code will not be submitted with the thesis 

% 
hlines = 14; 
inDir = '/Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Code (Github)/Icarus/inFiles/';
inDataFlame = '//Users/charlotteb/Desktop/Frontiers SIRS manuscript submission/Data/spectra/'; % Directory to spectra recorded with FLAME spectrometer at EOS

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
% [wcoh,wcs,f,coi] = wcoherence(s1i, s2i);
% Shift one against the other 
% a is 1 channel difference which is 0.07 nm 
% Each channel is 0.07 nm
s2icropA = s2i(2:end,1); wHgCropA = wHg(1:end-1,1); %
s1icropA = s1i(1:(end-1),1); 
% Sanity check 
figure
subplot(1,3,1)
plot(wHgCropA, s1icropA, '-k'); hold on; 
plot(wHgCropA, s2icropA, '-r'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
xlim([300 330]);

% b is 3 channel difference which is 0.4202
s2icropB = s2i(4:end,1); wHgCropB = wHg(1:end-3,1);
s1icropB = s1i(1:(end-3),1); 
% Sanity check 
subplot(1,3,2)
plot(wHgCropB, s1icropB, '-k'); hold on; 
plot(wHgCropB, s2icropB, '-r'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
xlim([300 330]);

% c is 6 channel difference which is 
s2icropC = s2i(7:end,1); wHgCropC = wHg(1:end-6,1);
s1icropC = s1i(1:(end-6),1); 
% Sanity check 
subplot(1,3,3)
plot(wHgCropC, s1icropC, '-k'); hold on; 
plot(wHgCropC, s2icropC, '-r'); 
xlabel('\lambda (nm)'); ylabel('Intensity (counts)'); 
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
xlim([300 330]);

figure
subplot(2,2,1)
% Magnitude-squared wavelet coherence 
% wcoh = wcoherence(x,y) returns the magnitude-squared wavelet coherence, which is a measure of the correlation between signals x and y in the time-frequency plane. 
% The inputs x and y must be equal length, 1-D, real-valued signals. The coherence is computed using the analytic Morlet wavelet.
% [wcoh, wcs] = wcoherence(x,y) returns the wavelet cross-spectrum of x and y. 
% You can use the phase of the wavelet cross-spectrum values to identify the relative lag between the input signals.
% [wcoh, wcs,f] = wcoherence(x,y,fs) uses the positive sampling frequency, fs, to compute the scale-to-frequency conversion, f. 
% The sampling frequency fs is in Hz.
% [wcoh,wcs,f,coi] = wcoherence(___) returns the cone of influence, coi, for the wavelet coherence in cycles per sample. If you specify the sampling frequency, fs, the cone of influence is in Hz.
% wcoherence(s1icrop, s2icrop, Fs, 'PhaseDisplayThreshold', 0.8); hold on; 
[~,wcsA,fwc,coiwc] = wcoherence(s1icropA, s2icropA, Fs); hold on; %Fs
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcsA(idxCoi) = NaN; 
% You can use the phase of the wavelet cross-spectrum values to identify the relative lag between the input signals.
magA = abs(wcsA);
thetaA = angle(wcsA); % The angles in theta are such that z = abs(z).*exp(i*theta). % phase angle 
% theta = angle(z) returns the phase angle in the interval [-π,π] for each element of a complex array z. The angles in theta are such that z = abs(z).*exp(i*theta).
% phase angle in the interval [-π,π] 
% degrees or radians that the waveform has shifted either left or right from the reference point.
% Plot phase as a function of frequency 
% phaseA = thetaA/pi; 
% pcolor(wHgCropA, fwc, phaseA); shading interp; colorbar; colormap(turbo(100)); 
% xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log'); ylim([min(fwc) max(fwc)]); % xlim([315.8 326.8]);
% cb = colorbar; 
% % cb.Limits = [0 3];
% ylabel(cb,'Phase/\pi', 'FontSize', 16, 'Rotation', -90); 
% cb.Label.Position(1) = 4;
% Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% Phase shift = 2 pi (T/T) where T is the common period 
degA = thetaA*180/pi;
pcolor(wHgCropA, fwc, degA); shading interp; colorbar; colormap(turbo(100)); 
xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); 
set(gca,'YScale','log'); ylim([min(fwc) max(fwc)]); xlim([min(wHgCropA) max(wHgCropA)]);
cb = colorbar;  WCS = abs(WCS);

% cb.Limits = [0 3];
ylabel(cb,'Phase (°)', 'FontSize', 12, 'Rotation', -90); 
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
% b 
[wcoh,wcsB,fwc,coiwc] = wcoherence(s1icropB, s2icropB, Fs); 
subplot(2,2,2)
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcsB(idxCoi) = NaN; 
magB = abs(wcsB);
thetaB = angle(wcsB); % The angles in theta are such that z = abs(z).*exp(i*theta). % phase angle 
degB = thetaB*180/pi;
pcolor(wHgCropB, fwc, degB); shading interp; colorbar; colormap(turbo(100));
xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log'); % xlim([315.8 326.8]); 
cb = colorbar; 
% cb.Limits = [0 6];
ylabel(cb,'Phase (°)', 'FontSize', 12, 'Rotation', -90);  
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
% c
subplot(2,2,3)
[wcoh,wcsC,fwc,coiwc] = wcoherence(s1icropC, s2icropC, Fs); 
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcsC(idxCoi) = NaN; 
magC = abs(wcsC);
thetaC = angle(wcsC); % The angles in theta are such that z = aCs(z).*exp(i*theta). % phase angle 
degC = thetaC*180/pi;
pcolor(wHgCropC, fwc, degC); shading interp; colorbar; colormap(turbo(100));
xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log'); % xlim([315.8 326.8]); 
cb = colorbar; 
% cC.Limits = [0 6];
ylabel(cb,'Phase (°)', 'FontSize', 12, 'Rotation', -90);  
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');

valA = nanmean(degA, 'all')

valB = nanmean(degB, 'all')

valC = nanmean(degC, 'all')

figure
plot([1,3,6], abs([valA, valB, valC]), '-ok')
xlabel('Shift applied (channels)'); ylabel('Mean phase (abs°)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

valA = nanmedian(degA, 'all')

valB = nanmedian(degB, 'all')

valC = nanmedian(degC, 'all')

figure
plot([1,3,6], abs([valA, valB, valC]), '-ok')
xlabel('Shift applied (channels)'); ylabel('Median phase (abs°)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');

% Adds window which will be extracted for figure 4 
maxF = 0.004;
minF = 0.001;
minW = 310;
maxW = 326.8; 
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

degA = degA(indxEFh:indxEFl,indxEWl:indxEWu); 
valA = nanmean(degA, 'all')

degB = degB(indxEFh:indxEFl,indxEWl:indxEWu); 
valB = nanmean(degB, 'all')

degC = degC(indxEFh:indxEFl,indxEWl:indxEWu); 
valC = nanmean(degC, 'all')

figure
plot([1,3,6], abs([valA, valB, valC]), '-ok')
xlabel('Shift applied (channels)'); ylabel('Mean phase (abs° in extraction window)');
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');






