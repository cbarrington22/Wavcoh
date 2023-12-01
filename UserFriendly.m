% This USER-FRIENDLY SCRIPT contains code to determine the minimum and differential magnitude-squared wavelet coherence (MSWC) between spectra recorded by scanning UV-instruments using 
% the approach outlined in “Barrington et al. (2023) Simplifying the retrieval of volcanic SO2 from UV spectra using wavelet coherence”. 
 
% Created by C. BARRINGTON, March 2023.

% DESCRIPTION of script 
% Sections of the code which are contained in the "Manuscript" folder at https://github.com/cbarrington22/Wavcoh have been re-written to provide a user-friendly option for those wishing 
% to apply similar analysis to their own UV spectra. Although this code is intended as a stand-alone script it may require adapations depending on the instrument configuration or user's needs. 
% As default, this code is set-up to analyse spectra recorded by the NOVAC (Network for Observation of Volcanic and Atmospheric Change) and determines the MSWC for all measurement spectra recorded within 
% a full instrument scan. It assumes the wavelength information contained in 'fnHg' (below) is correct. If that is not the case, but a differential SO2 absorption signal is observed in the MSWC, 
% the user should re-define the extraction window based on the differential SO2 absoprtion signal (which is expected to be centered at ~310 nm) (minWe, maxWe). Simillarly, the extraction window 
% is defined using the same frequency limits as indicted in the manuscript. In the case that the spectrometer resolution varies from that of instruments belonging to the NOVAC network, the frequency 
% limits (maxF and minF) should may need to be djusted in order to incorporate the equivalent spatial frequency range. This script differs slightly from the code contained in the "Manuscript" folder 
% at https://github.com/cbarrington22/Wavcoh but the analysis procedure is the same.

%% HOW TO USE THIS SCRIPT 
% This script is divided into three parts. 
% The first part defines the variables, many of which are defined as default, and must be run before either the second or third sections. 
% The second part of this script applies the approach outlined in Barrington et al. (2023) Simplifying the retrieval of volcanic SO2 from UV spectra using wavelet coherence, to define the 
% minimum and differential MSWC between spectra. It produces two figures, one with two subplots showing the minimum and differential MSWC respectively, and a second figure which displays the mean 
% phase angle of the wavelet cross spectrum (WCS) (both concerning only data within the extraction window). 
% The third part of this script applies the same approach but plots the MSWC and WCS over the entire wavelength and spatial frequency range. The purpose of this section is to allow users to verify the 
% suitability of the extraction window and futhur explore the differences in the MSWC and WCS which are observed.  

%% FORMAT OF SPECTRAL FILES 
% This script is set up to analyse spectral files which contain 3 columns, where: 
    % Column 1 - contains the spectral information and the name of the variable or information which is also recorded, either before or after the spectral information 
    % Column 2 - contains the '=' asigning the variable in column 1 to the data contained in column 3 information in column 3 to the variable contained in column 1 (and is redundant) 
    % Column 3 - contans data regarding the recorded spectrum e.g., the scan elevation angle  
% This is the configuration which enables .STD files recorded by the NOVAC to be read. Please amend lines 76 to 81 if a different format is required.  

%% PART 1: USER DEFINED VARIABLES 
%% MUST BE CHANGED
toAnalyse=''; % Defines the directory to the instrumnet scan
% to be analysed.
outDir=''; % Defines directory to save figures  
% Defines row indicies for intenisty data: 
firstChannel=4; % Enter the row number which corresponds to the FRIST LINE of spectral data 
lastChannel=2051; % Enter the row number which corresponds to the LAST LINE of spectral data 
nameInd=2091; % Enter the row number which corresponds to the name, ID or reference for the spectrum
% Note that this code is currently set up to analyse all spectra labelled as 'scan'. Please alter the switch statement starting at line 82 if a different ID should be used. 
angleInd=2072; % Enter the row number which corresponds to the elevation angle of the spectra. This code is set up to analyse spectra within -80 and 80 degrees from zenith (see 'limAng' below).
chanWav=''; % Directory to file containing wavelength information
fnHg=''; % Name of file containing wavelength information. This file should be .txt file with two columns with the same length as the measurement spectra
hlines = 0; % Enter the number of header lines to skip in 'fnHg' file (above)
%% MAY BE CHANGED
targrefAng=-60; % Scan elevation angle to be used as reference (note the closest angle may be used within the range -80 to 80 (limAng) and if 0 is selected, it is the measured spectra recorded 
% at zenith not the 'sky' spectrum typically used in the NOVAC algorithm. The 'sky' spectrum is not used in this script. (Default = -60).
limAng=80; % Enter the maximum elevation angle from zenith with which spectra should be analysed. Note that this value should always be positive. At scan elevation angles >80 degrees from zenith
% the radiative transfer effects (RTEs) are large, that changes in the MSWC at scan elevation angles beyond 80 degrees from zenith, may not be due to the differential absorption by SO2. (Default 
% = 80).
% Sets limits for extarction window. The assigned values depend on the wavelength range and resolution of the spectrometer. 
% Spatial frequency limits (cycles/wavelength) 
maxF=0.004;
minF=0.001;
% Wavelength range (nm)
minWe=310.0;
maxWe=326.8; 
% If the number of measurement spectra or the elevation angle range changes, the division of figures and subplots in part 3 will need changing. With the default configuration and spectra recorded 
% by the NOVAC, a total of 45 spectra will be analysed per instrument scan. To plot the MSWC and WCS for each spectrum, two figures are used for each (panel 1 and panel 2) each of which have space
% for 25 subplots. If the total number of spectra to be analysed exceeds 50, the definition of subplots in lines 216 and 236 (for MSWC) and 260 and 282 (for WCS) should be changed. 

%% PART 2: APPLIES APPROACH 
subFolders=dir(toAnalyse); % Returns data on subfolders to be analysed
subFoldersData=subFolders(~ismember({subFolders(:).name},{'.','..','.DS_Store'})); % Removes '.' and '..' and '.DS_Store' contained in subFolder varaible
spectraList=extractfield(subFoldersData,'name');spectraList=spectraList'; % Extracts file names to list all spectra to be analysed
sky=[]; % Creates variables used in loop (below) 
scan=[];
angleOut=[];

% EXTRACTS SPECTRAL INFORMATION from each file
for a=1:length(spectraList) % For all spectra in the instrument scan
    spectraName=spectraList(a); 
    spectraDirectory=fullfile(toAnalyse,spectraName); % Directory to access individual spectra within directory
    inSpectraI=char(spectraDirectory);
    [data1,~,data3]=textread(inSpectraI,'%s%s%s'); % Reads file
    spectra=data1(firstChannel:lastChannel,:); % Isolates intensity data 
    spectra=str2double(spectra); % Converts intenisty data to double 
    % Extracts iformation about recorded spectrum:
    name=data3(nameInd,:);name=string(name); % Name or ID
    aNgle=data3(angleInd,:); % Scan angle 
    switch name % Compiles spectral information from all measurement spectra
        case '"sky"'
        sky=[sky,spectra]; % Stores intenisty data from 'sky' spectrum in case the user wants to find the differential absorption signature of SO2 between the measurement spectra recorded at 
        % zenith and the 'sky' spectrum  
        case '"dark"' % Ignores dark 
        case '"scan"'
        scan=[scan,spectra];
        angleOut=[angleOut,aNgle]; % Compiles corresponding scan elevation angle 
        otherwise 
        warning('No matching case for line %d in the spectral file.',nameInd); % Warns user that the spectrum ID does not match one of the cases (check nameInd and lines 82 to 87)
    end
end 

% IDENTIFIES SPECTRUM TO BE USED AS THE REFERENCE 
scanAngle=str2double(angleOut); % NOVAC full instrument scans record 51 measurement spectra 
[~,idx]=min(abs(scanAngle-targrefAng)); % Reference angle index
refAng=scanAngle(idx); % Scan elevation angle of reference spectrum
refSpec=scan(:,idx); % Spectral information for reference 
% Find scan angles between range of interest (within defined limits)
[~,cSpec]=find((scanAngle<limAng)&scanAngle>(-limAng));
scanAngleAn=scanAngle(:,cSpec); % For NOVAC instruments, scanAnglesAn should contain angles -79 to 79 degrees from zenith if using 80 degrees from zenith as limAng
scanAn=scan(:,cSpec); % Isolates spectra to be analysed (these are the 45 useful scans to be analysed and includes the reference spectra as a sanity check) 

% DETERMINES LIMITS OF EXTRACTION WINDOW 
dHg = fullfile(chanWav,fnHg); inHg = char(dHg); % Loads Hg spectrum or wavelength information
if hlines >= 1  % Skips header lines if exist 
   [wHg, iHg] = textread(inHg,'%f %f','headerlines',hlines); % Isolates wavelength information 
else
   [wHg, iHg] = textread(inHg,'%f %f'); 
end  

% DETERMINES SAMPLING FREQUENCY (Fs)
chnls=length(spectra); % Determines number of spectrometer channels   
maxW=max(wHg); minW=min(wHg); % Determines minimum and maximum wavelength
rangeW=maxW-minW; % Determines wavelength range of spectrometer 
Fs=rangeW/chnls; % Fs is the sampling rate in wavelength, each data point (or channel) is equivalent to Fs wavelengths (approximately) - note this is only used to determine spatial frequency in cycles/lambda 

% DEFINES EXTRACTION WINDOW 
[wcoh,~,fwc,coiwc]=wcoherence(refSpec,scanAn(:,1),Fs); % Determines MSWC, frequency scales and cone of influence (COI) 
% Wavelength
% Finds channel (wavelength row index) which is closest to the upper and lower wavelength limits of the analysis window 
[~,indxEWl]=min(abs(wHg-minWe)); % Lower 
WlAw=wHg(indxEWl); 
[~,indxEWu]=min(abs(wHg-maxWe)); % Upper
WuAw=wHg(indxEWu); 
% Frequency
% Frequency levels in F are descending, so low row indices represent higher frequencies and higher row indices represent low frequencies 
[~,indxEFh]=min(abs(fwc-maxF)); % Highest 
F1Aw=fwc(indxEFh); 
[~,indxEFl]=min(abs(fwc-minF)); % Lowest 
F2Aw=fwc(indxEFl); 
% Finds indicies for COI 
szCoi=length(coiwc);szF=length(fwc);matF=repelem(fwc,1,szCoi);colCoi=coiwc';matCoi=repelem(coiwc,szF,1);idxCoi=matF<=matCoi; % Defines index for cone of influence (COI)

% ALLOCATES SPACE to variables used in analysis loop 
wcohEW=wcoh(indxEFh:indxEFl,indxEWl:indxEWu); % MSWC in extraction window 
wcohV=wcohEW(:); % Vectorises
wcohEx=NaN(length(wcohV),width(scanAn)); % Extracted MSWC
wcsEx=NaN(length(wcohV),width(scanAn)); % Extracted WCS
minMSWCe=NaN(1,width(scanAn)); % Minimum MSWC
sumMSWCe=NaN(1,width(scanAn)); % Sum of MSWC
paWCS=NaN(1,width(scanAn)); % Mean phase angle of WCS 

% RUN ANALYSIS to determine MSWC and WCS
for i=1:width(scanAn) % For all spectra considered 
    y=scanAn(:,i);
    [wcoh,wcs,fwc,~]=wcoherence(refSpec,y,Fs); % Determine MSWC and WCS between spectra and reference 
    wcoh(idxCoi)=NaN; % Removes data within COI
    wcs(idxCoi)=NaN; % Removes data within COI
    c=abs(wcoh(indxEFh:indxEFl,indxEWl:indxEWu)); % Extract MSWC within extraction window 
    wcohEx(:,i)=c(:); % Save 
    x=abs(wcs(indxEFh:indxEFl,indxEWl:indxEWu)); % Extract WCS within extraction window 
    wcsEx(:,i)=x(:); % Save 
end 

% FIGURES 
% PLOTS MSWC (Figs. 6, 7 and 9 in manuscript) 
minMSWC=1-(min(wcohEx)); % Finds minimum MSWC in extraction window 
bl=5.44e3; % Defines sum of MSWC in extraction window for two 'clear-sky' spectra (taken from Fig. 5 in manuscript). Note: This value is taken from MSWC between spectra with slightly different 
% wavelength range and Fs. This value is therefore only an approximation for the sum of the MSWC between two 'clear-sky' spectra. The purpose however is to enable a measure of 'differential' MSWC so the trend is intuitive. 
diffMSWC=(sum(wcohEx))-bl; % Finds differential MSWC in extraction window 
figure
% Minimum MSWC
subplot(1,2,1)
p=plot(scanAngleAn,minMSWC,'-ok'); hold on; % Plots minimum MSWC against elevtaion angle of scan
p.MarkerFaceColor='k';
xlabel('Scan elevation angle (° from zenith)');
ylabel('1 - minimum MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig=gca;Fig.FontSize=12;set(gcf,'color','w');
xl=xline(refAng,':');xl.LineWidth=2;
xlim([-(limAng+3) limAng+3])
% Differential MSWC
subplot(1,2,2)
p=plot(scanAngleAn,abs(diffMSWC),'-ok'); hold on; % Plots differential MSWC against scan elevation angle
p.MarkerFaceColor='k';
xlabel('Scan elevation angle (° from zenith)');
ylabel('Differential MSWC'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig=gca;Fig.FontSize = 12;set(gcf,'color','w');
xl=xline(refAng,':');xl.LineWidth=2;
xlim([-(limAng+3) limAng+3])
fnOut='MSWC.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 

% PLOTS MEAN PHASE ANGLE OF WCS (Figs. S8 and S10 in manuscript) 
thetaX=angle(wcsEx); % Returns the phase angles (in radians) of the WCS for each spectra analysed 
degX=thetaX*180/pi; % Converts units to degrees
meanPA=nanmean(degX); % Takes the mean 
stdPA=nanstd(degX); % Finds standard deviation 
figure
p=errorbar(scanAngleAn,meanPA,stdPA,stdPA,'vertical','-ok'); hold on; % Plots mean phase angle against scan elevation angle
p.MarkerFaceColor='g';
xlabel('Scan elevation angle (° from zenith)');
ylabel('Mean phase (°)'); 
set(gca,'XMinorTick','on','YMinorTick','on'); 
Fig=gca;Fig.FontSize=12;set(gcf,'color','w');
xl=xline(refAng,':'); xl.LineWidth=2;
xlim([-(limAng+3) limAng+3])
fnOut='Shift.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 

%% PART 3: MSWC AND WCS PANEL PLOTS 
% Defines display variables 
maxWC=0.9; % Maximum MSWC to be displayed 
Fh=maxF-minF; 
posE=[minWe,minF,maxWe-minWe,Fh]; % Dimensions for rectangle (to display extraction window)
% MSWC 
% Panel 1
figure
set(gcf,'units','points','position',[100,100,1500,1500])
for i=1:round(width(scanAn)/2)+2 % For first half of spectra 
    y=scanAn(:,i);
    [wcoh,~,fwc,~]=wcoherence(refSpec,y,Fs); % Determine MSWC and WCS between spectra and reference 
    wcoh(idxCoi)=NaN; % Removes data within COI
    subplot(5,5,i)
    pcolor(wHg,fwc,abs(wcoh));shading interp;colormap(hot(50));set(gca,'YScale','log'); hold on; % Plots MSWC
    rectangle('Position', posE,'LineWidth', 1.5,'LineStyle', '-', 'EdgeColor', 'k'); % Plots extraction window 
    box on
    Ang=string(scanAngleAn(:,i));
    title(sprintf('%s°',Ang));
    Fig=gca;Fig.FontSize=7;box on 
    xlabel('\lambda (nm)')
    ylabel('Spatial frequency (cycles/\lambda')
    Fig=gca;Fig.FontSize=10;set(gcf,'color','w');
end
cb=colorbar;cb.Limits=[0 maxWC];fnOut='MSWC_panel1.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 
% Panel 2 
figure
set(gcf,'units','points','position',[100,100,1500,1500])
run=1; 
for i=i+1:width(scanAn) % For remaining spectra  
    y=scanAn(:,i);
    [wcoh,~,fwc,~]=wcoherence(refSpec,y,Fs); % Determine MSWC and WCS between spectra and reference 
    wcoh(idxCoi)=NaN; % Removes data within COI
    subplot(5,5,run)
    pcolor(wHg,fwc,abs(wcoh));shading interp;colormap(hot(50));set(gca,'YScale','log'); hold on; % Plots MSWC
    rectangle('Position', posE,'LineWidth', 1.5,'LineStyle', '-', 'EdgeColor', 'k'); % Plots extraction window 
    box on
    Ang=string(scanAngleAn(:,i));
    title(sprintf('%s°',Ang));
    Fig=gca;Fig.FontSize=7;box on 
    xlabel('\lambda (nm)')
    ylabel('Spatial frequency (cycles/\lambda')
    Fig=gca;Fig.FontSize=10;set(gcf,'color','w');
    run=run+1;
end
cb=colorbar;cb.Limits=[0 maxWC];fnOut='MSWC_panel2.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 

% WCS 
% Panel 1
figure
set(gcf,'units','points','position',[100,100,1500,1500])
for i=1:round(width(scanAn)/2)+2 % For first half of spectra 
    y=scanAn(:,i);
    [~,wcs,fwc,~]=wcoherence(refSpec,y,Fs); % DetermineWCS between spectra and reference 
    wcs(idxCoi)=NaN; % Removes data within COI
    thetaA = angle(wcs); 
    degA = thetaA*180/pi;
    subplot(5,5,i)
    pcolor(wHg,fwc,degA);shading interp;colormap(turbo(50));set(gca,'YScale','log'); hold on; % Plots WCS
    rectangle('Position', posE,'LineWidth', 1.5,'LineStyle', '-', 'EdgeColor', 'k'); % Plots extraction window
    box on
    Ang=string(scanAngleAn(:,i));
    title(sprintf('%s°',Ang));
    Fig=gca;Fig.FontSize=7;box on 
    xlabel('\lambda (nm)')
    ylabel('Spatial frequency (cycles/\lambda')
    Fig=gca;Fig.FontSize=10;set(gcf,'color','w');
end
colorbar;fnOut='WCS_panel1.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 
figure
set(gcf,'units','points','position',[100,100,1500,1500])
run=1; 
for i=i+1:width(scanAn) % For remaining spectra  
    y=scanAn(:,i);
    [wcoh,wcs,fwc,~]=wcoherence(refSpec,y,Fs); % Determine MSWC and WCS between spectra and reference 
    wcs(idxCoi)=NaN; % Removes data within COI
    thetaA = angle(wcs); 
    degA = thetaA*180/pi;
    subplot(5,5,run)
    pcolor(wHg,fwc,degA);shading interp;colormap(turbo(50));set(gca,'YScale','log'); hold on; % Plots WCS
    rectangle('Position', posE,'LineWidth', 1.5,'LineStyle', '-', 'EdgeColor', 'k'); % Plots extraction window
    box on
    Ang=string(scanAngleAn(:,i));
    title(sprintf('%s°',Ang));
    Fig=gca;Fig.FontSize=7;box on 
    xlabel('\lambda (nm)')
    ylabel('Spatial frequency (cycles/\lambda')
    Fig=gca;Fig.FontSize=10;set(gcf,'color','w');
    run=run+1;
end
cb=colorbar;fnOut='WCS_panel2.png';fname=fullfile(outDir,fnOut);saveas(gcf,fname) % Saves figure accoridng to output directory 

% END OF SCRIPT    
