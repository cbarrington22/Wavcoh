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

% PLOTS RESAMPLED AND CONVOLUTED ABSORPTION CROSS SECTIONS AND SOLAR REFERENCE 
colsolar = [0 0 0];
figure 
dp = plot(wHg, cSolar, '-'); hold on; dp.LineWidth = 1; dp.Color = colsolar;
ylim([0 max(cSolar)]);
ylabel('Intensity (counts)'); xlim([min(wHg) max(wHg)]); xlabel('\lambda (nm)');
% dp = xline(lAw, ':k'); dp.LineWidth = 1; % Places vertical line at lower wavelength range of analysis window 
% dp = xline(uAw, ':k'); dp.LineWidth = 1; % Places vertical line at upper wavelength range of analysis window 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
 
% DETERMINES SAMPLING FREQUENCY (Fs)
% [...,F] = cwt(...,Fs) specifies the sampling frequency, Fs, in hertz as a positive scalar and returns the scale-to-frequency conversions in hertz, F. 
% If the sampling frequency is not specified, cwt returns F in cycles/sample (i.e. channel). 
% If the sampling frequency is specified (in wavelength), cwt returns F in cycles/wavelength. 
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

% Figure 6
% a
% Intensity of two measured sky spectra at zenith 
inSpectraS = '90_0_FLMS195681__0__10-13-49-165.txt'; % Sky zenith 
inSpectraSlea = '30_0_FLMS195681__0__11-26-46-904.txt'; % Sky large elevation angle 
sS = fullfile(inDataFlame, inSpectraS); sS = char(sS); sSlev = fullfile(inDataFlame, inSpectraSlea); sSlev = char(sSlev);
if hlines >= 1
   [~, sSi] = textread(sS,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
   [~, sSlevi] = textread(sSlev,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [~, sSi] = textread(sS,'%f %f'); 
   [~, sSlevi] = textread(sSlev,'%f %f'); 
end  
% Magnitude-squared wavelet coherence 
% Figure 6 
figure  
% a 
maxF = 0.004; % Copied from Fig3to6.m 
minF = 0.001;
minW = 310;
maxW = 326.8; 
subplot(1,2,1) 
[wcoh, ~, fwc, coiwc] = wcoherence(sSi, sSlevi, Fs);
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI
wcoh(idxCoi) = NaN; 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([280 max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
set(gca,'XMinorTick','on','YMinorTick','on'); 
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 16, 'Rotation', -90); 
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% Adds wavelength ranges of analysis window
Fh = maxF-minF; 
posE = [minW, minF, maxW-minW, Fh];
rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );
% b 
subplot(1,2,2) 
wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
xlim([minW maxW]);xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
set(gca,'XMinorTick','on','YMinorTick','on'); ylim([minF maxF]);
cb = colorbar; 
cb.Limits = [0 1];
ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 16, 'Rotation', -90); 
cb.Label.Position(1) = 4;
Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
val = sum(wcohEx,'all'); 


% % Figure 7 
% % IMPORTANT - wavelength ranges in this figure relate to those used in DOAS (not in extracting SO2 signal) 
% figure 
% % a
% subplot(2,2,[1 2])
% [wcoh, ~, fwc, coiwc] = wcoherence(cSolar, sSi, Fs);
% wcoh(idxCoi) = NaN; 
% pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
% xlim([280 max(wHg)]); xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
% cb = colorbar; 
% cb.Limits = [0 1];
% ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 16, 'Rotation', -90); 
% cb.Label.Position(1) = 4;
% Fig = gca; Fig.FontSize = 15; set(gcf,'color','w');
% val = sum(wcohEx,'all'); 
% % Primary and secondary wavelength ranges 
% lAw1 = 314.8;
% uAw1 = 326.8; 
% lAw2 = 326.8; 
% uAw2 = 353.5; 
% % Adds wavelength ranges of analysis window
% Fh = max(fwc)-min(fwc); 
% pos1wc = [lAw1, min(fwc), uAw1-lAw1, Fh];
% pos2wc = [lAw2, min(fwc), uAw2-lAw2, Fh];
% rectangle('Position', pos1wc,'LineWidth', 1,'LineStyle', '-', 'EdgeColor', [0, 0, 0, 0.6] );
% rectangle('Position', pos2wc,'LineWidth', 1,'LineStyle', '--', 'EdgeColor', [0, 0, 0] ); 
% % 
% subplot(2,2,3)
% % b % This time we don't limit the frequency range 
% % Adds window which will be extracted for figure 4 
% maxF = max(fwc);
% minF = min(fwc);
% % maxF = 0.004;
% % minF = 0.001;
% % minW = 310;
% % maxW = 326.8; 
% % DEFINES ANALYSIS WINDOW
% % WAVELENGTH 
% % Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
% [~, indxEWl] = min(abs(wHg - lAw1)); % Lower 
% WlAw = wHg(indxEWl); 
% [~, indxEWu] = min(abs(wHg - uAw1)); % Upper
% WuAw = wHg(indxEWu); 
% % FREQUENCY
% % Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
% [~, indxEFh] = min(abs(fwc - maxF)); % Highest 
% F1Aw = fwc(indxEFh); 
% [~, indxEFl] = min(abs(fwc - minF)); % Lowest 
% F2Aw = fwc(indxEFl); 
% % fwcEx = fwc(indxEFh:indxEFl, :); 
% % wHgEx = wHg(indxEWl:indxEWu, :); 
% % minMSWC = zeros(4, 1); 
% % minMSWC(1,1) = min(wcohEx,[],'all'); % Returns the smallest element
% % % Adds wavelength ranges of analysis window
% % Fh = maxF-minF; 
% % posE = [lAw1, minF, uAw1-lAw1, Fh];
% % rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );
% % posE = [lAw2, minF, uAw2-lAw2, Fh];
% % rectangle('Position', posE,'LineWidth', 0.5,'LineStyle', '-', 'EdgeColor', [1, 1, 1, 1] );
% % wcohEx = abs(wcoh(:,indxEWl:indxEWu)); 
% pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
% xlim([lAw1 uAw1]);xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
% set(gca,'XMinorTick','on','YMinorTick','on'); % ylim([minF maxF]);
% cb = colorbar; 
% cb.Limits = [0 1];
% ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 12, 'Rotation', -90); 
% cb.Label.Position(1) = 4;
% Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
% wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
% val1 = nanmean(wcohEx, 'all');
% % c % This time we don't limit the frequency range 
% % Finds channel (wavelength row index) which is closest to upper and lower wavelength limits of the analysis window 
% [~, indxEWl] = min(abs(wHg - lAw2)); % Lower 
% WlAw = wHg(indxEWl); 
% [~, indxEWu] = min(abs(wHg - uAw2)); % Upper
% WuAw = wHg(indxEWu); 
% % FREQUENCY
% % Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
% [~, indxEFh] = min(abs(fwc - maxF)); % Highest 
% F1Aw = fwc(indxEFh); 
% [~, indxEFl] = min(abs(fwc - minF)); % Lowest 
% F2Aw = fwc(indxEFl); 
% subplot(2,2,4)
% pcolor(wHg, fwc, abs(wcoh)); shading interp; colorbar; colormap(turbo(100)); 
% xlim([lAw2 uAw2]);xlabel('\lambda (nm)'); ylabel('Spatial frequency (cycles/\lambda)'); set(gca,'YScale','log');
% set(gca,'XMinorTick','on','YMinorTick','on'); % ylim([minF maxF]);
% cb = colorbar; 
% cb.Limits = [0 1];
% ylabel(cb,'Magnitude-squared wavelet coherence', 'FontSize', 12, 'Rotation', -90); 
% cb.Label.Position(1) = 4;
% Fig = gca; Fig.FontSize = 11; set(gcf,'color','w');
% wcohEx = abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); 
% val2 = nanmean(wcohEx, 'all');
