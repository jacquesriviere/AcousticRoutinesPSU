function [MaxInter,TimeShift,RmsAmp,Amp,RmsAmpRef,AmpRef,fullWFref,LocalAcTime] = ...
    ProcessAc_Tomo(WF_path,AcSettingsfile,SyncFile,idxBeg,idxEnd, ...
              reference,NtoStack,threshold,Offset,displayoptions,TimeRange)

% ProcessAc_Tomo cross-correlate all WFs with a reference WF to
% extract a change in time of flight. It also computes the amplitude of
% each waveform (peak-peak amplitude) and the root-mean-square amplitude of
% each WF.
% Waveforms can be stacked to improve SNR by choosing 'NtoStack' larger
% than 1.
% When choosing 'absref', all waveforms within the time range are stacked
% to create a reference waveform. When choosing 'relref', each stacked WF
% is cross-correlated with the previous stacked WF, unless the previous WF
% is noisy (i.e. maxinter below a certain threshold). When choosing
% 'mixref', the first Ntostack WFs are stacked and used as a reference. The
% reference is changed only when the maxinter becomes lower than a certain
% threshold.

% INPUTS
% WF_path is the path where the acoustic data can be loaded
% AcSettingsfile is the path where pXXXX.mat (containing acoustic settings) can be
% loaded

% SyncFile is the result of the AcSync routine, providing the acTime vector
% and the adjusted sampling time ts returned by SyncAcData function

% idxBeg is the idx corresponding to the beginning of the waveform to be
% analyzed
% idxEnd is the idx corresponding to the end of the waveform to be
% analyzed

% reference is set to 'absref', 'relref' or 'mixref'. 

% threshold is used when using 'relref' or 'mixref'. For 'relref', the
% current WF won't be used as a reference for the next WF if maxinter is
% below the threshold (i.e. noisy waveform). For 'mixref', the current WF
% will become a reference if maxinter is below threshold.

% TimeShift and RmsAmp values are kept only when MaxInter is above this
% threshold. threshold is comprised between -1 (no threshold) and 1 (max
% threshold). When 'absref' is used, threshold is ignored.

% Offset is a number allowing one to offset waveforms
% corresponding to different channels when displaying waveforms

% displayoptions: if 0, all waveforms will be displayed. if 1, only 1 over
% 100 waveforms will be displayed

% TimeRange: Time range that will be analyzed. If not defined, the whole
% run will be analyzed

% OUTPUTS
% MaxInter is the maximum of intercorrelation vector

% TimeShift is the delay found between each WF with its reference WF
% (either absolute or relative reference). When the reference is relative,
% cumsum can be used to sum up the cumulative time differences.

% RmsAmp contains the rms amplitude of each WF

% RmsAmpRef is the rms amplitude of the reference WF 
% fullWFref is the reference WF.
% LocalAcTime is the new acoustic time vector in the considered TimeRange provided. 
% It accounts for stacking and the number of transmitters used during the
% analysis. The size of LocalActime is [number of stacked waveforms per transmitter; number of transmitters]

format short

if strcmp(reference,'absref')
    threshold = -1;
    warning('The threshold is set to -1 (i.e., no threshold) when using ''absref''.')    
end

% acoustic parameters
acSettings = load(AcSettingsfile);                          % load acoustic settings
numSFpfile = acSettings.numFrames/2;                        % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;                           % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;                     % number of WF per file and per channel
numCHR = length(acSettings.channels2save);                  % number of receiving channels
numCHT = length(acSettings.channels2transmit);              % number of transmitting channels
WFlength = acSettings.Nsamples;                             % waveform length
clear acSettings

load(SyncFile); % load sync data: acTime, ts, acPeriod, totalnumberoffiles

%% find the files and WFs within the files where analysis begins and ends
if nargin < 11 % if TimeRange is not specified, the whole run is analyzed.        
    filenumber1 = 1; % first file
    filenumber2 = totalnumberoffiles; % last file
    idxWFwithinfile1 = 1; % first WF to be analyzed within the first file
    idxWFwithinfile2 = numWFpfilepCH; % last WF to be analyzed within the last file
    idxacTime1 = 1;
    idxacTime2 = length(acTime);    
else    
    [filenumber1,idxWFwithinfile1,idxacTime1,acTime1] = findidxs(acTime,TimeRange(1),1,numCHT,numWFpfilepCH);
    [filenumber2,idxWFwithinfile2,idxacTime2,acTime2] = findidxs(acTime,TimeRange(2),numCHT,numCHT,numWFpfilepCH);     
end
display(['From WF ' num2str(idxWFwithinfile1) ' of file ' num2str(filenumber1) ' to WF ' num2str(idxWFwithinfile2) ' of file ' num2str(filenumber2) '.']);

acN = idxacTime2 - idxacTime1 + 1; % total number of WF per receiver
acN = floor(acN/NtoStack/numCHT); % total number of stacked WF per receiver and per transmitter. "floor" because NtoStack is not necessarily a multiple of the total number of waveforms

% define LocalAcTime to account for stacking, the number of transmitters
% used and the time range where the analysis is conducted

acPeriod_new = NtoStack*acPeriod*numCHT;
LocalAcTime = (0:acN-1)*acPeriod_new/1e6; %acPeriod is in microsec

acTime_newshifted = zeros(length(LocalAcTime),numCHT); % matrix containing one time vector per transmitter
for chnumt = 1:numCHT % shift by the average time of all stacked WFs
    acTime_newshifted(:,chnumt) = LocalAcTime + mean([acTime(idxacTime1+chnumt-1) acTime(idxacTime1+chnumt-1+(NtoStack-1)*numCHT)]);           
end

LocalAcTime = acTime_newshifted;

% time vector for each waveform
timeWF = (0:WFlength-1)*ts;

fullWFref = zeros(WFlength,numCHR,numCHT);
RmsAmpRef = zeros(numCHR,numCHT);
AmpRef = zeros(numCHR,numCHT);

% filter to be used if noisy waveforms (adjust order and frequencies as needed)
filterparam = fir1(256,[0.25 2]*ts*2); % pass band filter 0.25MHz 2MHz (ts is in microsec)

%% build a reference waveform

% if 'absref' is chosen, all waveforms in the time range are stacked
% to build a template

% if 'relref' or 'mixref' is chosen, this reference WF will be used only once to
% be compared with the next one

kk = 0; % from 0 to NtoStack*numCHT - 1 for relref and mixref or from 0 to acN*numCHT - 1 fr absref 
ii = filenumber1; % file number
jj = idxWFwithinfile1; % from 1 to numWFpfilepCH 
chnumt = 1; % transmitter index

if strcmp(reference,'relref')||strcmp(reference,'mixref')
    upperlimit = NtoStack*numCHT;
elseif strcmp(reference,'absref')
    upperlimit = acN*NtoStack*numCHT;   
end
while kk < upperlimit % number of Reference WFs    
    if jj == 1 || (ii == filenumber1 && jj == idxWFwithinfile1) % open new file if jj = 1 or if it's the first file
        ACdata = LoadAcFile(WF_path,ii,numCHR,numSFpfile);        
    end        
    fullWFref(:,:,chnumt) = fullWFref(:,:,chnumt) + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % stack WFs
    if chnumt < numCHT % chnumt runs from 1 to numCHT
        chnumt = chnumt + 1;
    else
        chnumt = 1;
    end   
    if jj < numWFpfilepCH  % stay within the same file for the next run
        jj = jj + 1;
    else                   % use next file for the next run
        jj = 1;ii = ii + 1;
    end
    kk = kk + 1;
end

if strcmp(reference,'relref')||strcmp(reference,'mixref')
    fullWFref = fullWFref/NtoStack; % average
elseif strcmp(reference,'absref')
    fullWFref = fullWFref/acN/NtoStack; % average
end

% figure(765);plot(fullWFref(:,1,1));hold on; % uncomment to display the effect of filtering
% fullWFref = filtfilt(filterparam,1,fullWFref); % filtering
% plot(fullWFref(:,1,1),'k');hold off;pause % uncomment to display the effect of filtering

WFref = fullWFref(idxBeg:idxEnd,:,:); % part of the WF to be analyzed

for chnumr = 1:numCHR
    RmsAmpRef(chnumr,:) = rms(WFref(:,chnumr,:)); % RmsAmp of the reference waveform
    AmpRef(chnumr,:) = max(WFref(:,chnumr,:))-min(WFref(:,chnumr,:)); % RmsAmp of the reference waveform
end

figure
idxoffset = 1;
for chnumt = 1:numCHT
    for chnumr = 1:numCHR
        %     subplot(numCHR,1,chnumr);
        plot(timeWF,fullWFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'r');hold on;grid on
        plot(timeWF(idxBeg:idxEnd),WFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'g');        
        legend('FullWF','WFtoAnalyze');
        xlabel('Time ($\mu s$)','Interpreter','Latex');
%         ylabel('Amplitude (bits)','Interpreter','Latex');
        
        dcmObj = datacursormode;
        set(dcmObj,'UpdateFcn',@GoodCursor);
        set(gca,'FontSize',16);
        set(gca,'YTickLabel',[])
        drawnow
        idxoffset = idxoffset + 1;
    end
end
hold off

clear ACdata

% return % uncomment here to look at the reference waveform

%% Compute changes in time of flight, max of intercorrelation and RmsAmp

MaxInter = zeros(acN,numCHR,numCHT);
TimeShift = zeros(acN,numCHR,numCHT);
RmsAmp = zeros(acN,numCHR,numCHT);
Amp = zeros(acN,numCHR,numCHT);

ii = filenumber1; % file number
jj = idxWFwithinfile1; % from 1 to numWFpfilepCH

% adjust depending on what you want to display
h1 = zeros(numCHR,numCHT);
h2 = zeros(numCHR,numCHT);
h3 = zeros(numCHR,numCHT);
h4 = zeros(numCHR,numCHT);

figure
set(gca,'YTickLabel',[]);
set(gca,'Ylim',[-numCHT*numCHR*Offset Offset]);
set(gca,'NextPlot','replacechildren');

for hh = 1:acN % from 1 to the total number of stacked waveforms
    
    fullWF = zeros(WFlength,numCHR,numCHT);    
    % stack WFs
    chnumt = 1; % transmitter index
    for kk = 1:NtoStack*numCHT  
        if (jj == 1) || (ii == filenumber1 && jj == idxWFwithinfile1) % open new file if jj = 1 or if it's the first file
            fprintf(['File number ' num2str(ii) '.\n']) % display file number            
            ACdata = LoadAcFile(WF_path,ii,numCHR,numSFpfile); 
        end
        fullWF(:,:,chnumt) = fullWF(:,:,chnumt) + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % read data
        if chnumt < numCHT % chnumt runs from 1 to numCHT
            chnumt = chnumt + 1;
        else
            chnumt = 1;
        end        
        
        if jj < numWFpfilepCH   % stay within the same file for the next run
            jj = jj + 1;
        else                    % use next file for the next run
            jj = 1;ii = ii + 1;
        end
    end
    fullWF = fullWF/NtoStack; % stacked WF  
    
%     figure(765);plot(fullWF(:,1,1));hold on; % uncomment to display the effect of filtering    
%     fullWF = filtfilt(filterparam,1,fullWF); % filtering    
%     plot(fullWF(:,1,1),'k');hold off;pause % uncomment to display the effect of filtering
    
    WF = fullWF(idxBeg:idxEnd,:,:); % WF is only the part to be analyzed    
    
    % cross-correlate (time delay)
    for chnumt = 1:numCHT
        for chnumr = 1:numCHR                        
            corr_signals = xcorr(WFref(:,chnumr,chnumt),WF(:,chnumr,chnumt),'coeff');
            [MaxInter(hh,chnumr,chnumt),TimeShift(hh,chnumr,chnumt)] = delay(corr_signals,ts);                                                                                        
        end
    end
    % amplitudes    
    RmsAmp(hh,:,:) = rms(WF,1);               % RmsAmp of the waveform
    Amp(hh,:,:) = max(WF,[],1)-min(WF,[],1);  % Max Amp of the waveform            
    
    if hh/100==round(hh/100) % display Max intercorrelation every 100 stacked waveforms
        MaxInterCorrelation = squeeze(MaxInter(hh,:,:));
        display(MaxInterCorrelation)
    end
    
    % display waveforms        
    if (displayoptions == 0) || (displayoptions == 1 && (hh-1)/100==round((hh-1)/100))
        idxoffset = 1;        
        for chnumt = 1:numCHT
            for chnumr = 1:numCHR                 
                h1(chnumr,chnumt) = plot(timeWF,fullWFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'r');hold on                
                h2(chnumr,chnumt) = plot(timeWF(idxBeg:idxEnd),WFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'g'); %hold on;                
                h3(chnumr,chnumt) = plot(timeWF,fullWF(:,chnumr,chnumt)-Offset*(idxoffset-1),'b');
                h4(chnumr,chnumt) = plot(timeWF(idxBeg:idxEnd),WF(:,chnumr,chnumt)-Offset*(idxoffset-1),'k'); %hold on
                drawnow
                idxoffset = idxoffset + 1; %pause         
            end
        end
        cla(gca); % clear former WFs without removing axis properties
    end   

    % choose reference WF for the next run when relref is used.
    for chnumt = 1:numCHT
        for chnumr = 1:numCHR
            if strcmp(reference,'relref') && MaxInter(hh,chnumr,chnumt) >= threshold                
                % current WF becomes reference WF for the next run if max
                % of intercorrelation is above threshold
                fullWFref(:,chnumr,chnumt) = fullWF(:,chnumr,chnumt);
                WFref(:,chnumr,chnumt) = WF(:,chnumr,chnumt);                
            elseif strcmp(reference,'relref') && MaxInter(hh,chnumr,chnumt) < threshold                
                % MaxInter is lower than threshold.
                % result of cross-correlation is replaced by NaN and we
                % don't keep the current WF as a reference for the next run
                TimeShift(hh,chnumr,chnumt) = NaN;
                RmsAmp(hh,chnumr,chnumt) = NaN;                
            end
            % if mixref and last 5 waveforms have MaxInter below threshold 
            if hh >= 5 
                if strcmp(reference,'mixref') && sum(MaxInter(hh-4:hh,chnumr,chnumt) < threshold) == 5  % last 5 waveforms are below threshold
%                 if strcmp(reference,'mixref') && sum(MaxInter(hh,chnumr,chnumt) < threshold) == 1 % last waveform are below threshold
                    % current WF becomes reference WF for the next runs
                    fullWFref(:,chnumr,chnumt) = fullWF(:,chnumr,chnumt);
                    WFref(:,chnumr,chnumt) = WF(:,chnumr,chnumt);
                    
                    % keep track of when the ref has changed (use structure because the number of ref changes are not equal for all combinations of T-R)
                    chrefname = ['R' num2str(chnumr) 'T' num2str(chnumt)];
                    try
                        changeref.(chrefname) = [changeref.(chrefname) hh]; % keep track of when the ref has changed
                    catch
                        changeref.(chrefname) = hh;
                    end
                end
            end            
        end
    end
end

% take the opposite such that a positive TimeShift means later arrival.
TimeShift = -TimeShift;

% sum up all time differences when 'relref' is used
if strcmp(reference,'relref')
    idxNaN = isnan(TimeShift); % find NaNs
    TimeShift(idxNaN) = 0; % 0 instead of NaN to be able to use cumsum below
    
    % sum up all the delays and put back NaN into the final vector
    TimeShift = cumsum(TimeShift);
    TimeShift(idxNaN) = NaN;
end
% adjust TimeShift based on when the reference has changed (mixref options)
if strcmp(reference,'mixref') && exist('changeref','var')
    for chnumt = 1:numCHT
        for chnumr = 1:numCHR
            chrefname = ['R' num2str(chnumr) 'T' num2str(chnumt)];                        
            if isfield(changeref,(chrefname)) % test if the ref was changed for this pair of T & R            
                kk = 1;
                for gh = changeref.(chrefname)
                    if gh == changeref.(chrefname)(end) % if last index, correct until the end
                        upbound = length(TimeShift(:,chnumr,chnumt));
                    else
                        upbound = changeref.(chrefname)(kk+1); % correct until next reference
                    end
                    TimeShift(gh+1:upbound,chnumr,chnumt) = TimeShift(gh+1:upbound,chnumr,chnumt) + TimeShift(gh,chnumr,chnumt);
                    kk = kk + 1;
                end
            end
        end
    end
end

end

