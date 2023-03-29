% wavelet coherence 
% This code is used to produce the figures included in 'wavelet coherence' section of Chapter 2 
% This code will not be submitted with the thesis 

% Masaya scans 
% Take one full instrument scan 
% Select the scan with 60 degree elevation angle as the reference 
% Wavelet coherence of each subsequent scan (to 60 degrees) 
% Include Hg spectrum to obtain wavelength information 
% We assume that spectra are calibrated according to wavelength information of Hg spectrum 
% Extract analysis window accord to previously used window 
% What is the wavelet coherence threshold identified by the zneith-60 degree elevation angle test 
% That's the baseline - we only consider wavelet coherence which crosses that threshold total wc 5.4215e+03 (Fig. 6 ) 
% Plot a pannel of subplots showing SO2 absorption signature 
% Extract the total wavelet coherence of the extraction window 
% Save total wavelet coherence and elevtaion angle of each scan 

inDir = '/Users/charlotteb/Documents/Frontiers SIRS manuscript submission/Data/D2J2375Florian/Example_Florian/'; % Directory to Hg spectrum for spectrometer used to recorded spectra in inDir
toAnalyse = '/Users/charlotteb/Documents/Frontiers SIRS manuscript submission/Data/D2J2375Florian/D2J2375/D2J2375_140319_1417_0/'; % Directory to full full instrument scan

% D2J2375_140319_1941_0/'; % Directory to full full instrument scan
% D2J2375_140319_1417_0

hlines = 0; 
% Loads Hg spectrum to obtain wavelet information consistent with that used in linear model  
fnHg = 'HG_D2J2375_0.txt'; % Hg spectrum 
dHg = fullfile(inDir, fnHg); inHg = char(dHg); % Hg spectrum 
if hlines >= 1
   [wHg, iHg] = textread(inHg,'%f %f', 'headerlines', hlines); % Skips header lines if exist 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  

subFolders=dir(toAnalyse); % Returns data on subfolders to be analysed
subFoldersData=subFolders(~ismember({subFolders(:).name},{'.','..','.DS_Store'})); % Removes '.' and '..' and '.DS_Store' contained in subFolder varaible
spectraList=extractfield(subFoldersData,'name');spectraList=spectraList'; % Extracts file names to list all spectra to be analysed
dark=[];
sky=[];
scan=[];
angleOut=[];

for a=1:length(spectraList) 
    spectraName=spectraList(a); 
    spectraDirectory=fullfile(toAnalyse,spectraName); % Directory to access individual spectra within directory
    inSpectraI=char(spectraDirectory);
    [data]=textread(inSpectraI,'%s');
    name=data(2150,:); % For some reason this is different % 2150
    name=string(name);
    angle=data(2093,:); 
    spectra=data(4:2051,:); 
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
chnls = length(wHg); % Number of spectrometer channels  
minW = min(wHg); % Starting wavelength 
maxW = max(wHg); % End wavelength 
rangeW = maxW - minW; % Determines wavelength range of spectrometer 
Fs = rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths 

% convert angle to number 
scanAngle = str2double(angleOut);

% Reference angle 
targrefAng = 0; 
% [~,idx] = min(abs(scanAngle - targrefAng));
% refAng = scanAngle(idx);
% refAng = scanAngle(9);
idx = 26 ; 
refSpec = scan(:, idx); % Reference scan 

% [~,cSpec] = find((scanAngle<refAng) & scanAngle>(-60)); 
[~,cSpec] = find(scanAngle>(-(targrefAng)) & scanAngle<(targrefAng+1)); 
% [~,cSpec] = size(scan);

[~, ~, fwc, ~] = wcoherence(refSpec, scan(:, idx+1), Fs);

run = 1;
maxF = max(fwc); % 0.001; % 0.005; % max(fwc) % 0.005; 
minF = min(fwc); % 0.001); % min(fwc) % 0.001;
minW = 300; % 310;
maxW = 353; 
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

[wcoh, ~, fwc, coiwc] = wcoherence(refSpec, scan(:, idx+1), Fs);
szCoi = length(coiwc); szF = length(fwc); matF = repelem(fwc, 1, szCoi); colCoi = coiwc'; matCoi = repelem(coiwc, szF, 1); idxCoi = matF <= matCoi; % Removes COI

% For plotting/saving 
temp = wcoh(indxEFh:indxEFl,indxEWl:indxEWu); 
temp = temp(:);
wcohExOut = zeros(2, length(cSpec)); 
maxWC = 1 % 0.8; % 0.9871; % from figure 6 
out = []; 

figure 
for i = 1:width(scan)
    y = scan(:, i);
    [wcoh, ~, ~, ~] = wcoherence(refSpec, y, Fs);
    wcoh(idxCoi) = NaN; 
    % subplot(7, 5, run)
    subplot(9, 6, run)
    pcolor(wHg, fwc, abs(wcoh)); shading interp; colormap(hot(100)); set(gca, 'YScale', 'log')

    xlim([minW maxW]); ylim([minF maxF]); 
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

% Figure 5 
% Differencial coherence 
% masaya1 = out;
masaya2 = out;

figure % All in window 
% x = scanAngle(:,cSpec); 
% plot(x, 1-masaya1, '-ok'); hold on; % Plots minimum wavelet coherend against scan angle 
plot(scanAngle, 1-masaya2, '-ob'); hold on; % Plots minimum wavelet coherend against scan angle 
xlabel('Scan angle (°)'); ylabel('1 - minimum MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig = gca; Fig.FontSize = 14; set(gcf,'color','w');
legend('Scan at 19:41 UTC', 'Scan at 14:17 UTC')
% 









