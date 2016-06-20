function [MaxInter,MinInter,TimeShift,RmsAmp,Amp,TOF_0,RmsAmpRef,AmpRef,fullWFref,acTime_new,integerShift] = ...
    ProcessAc_Tomo(Ac_path,run_ac_path,ts,totalnumberoffiles,idxBeg,idxEnd,idx_TOF_0, ...
              stepoptions,reference,NtoStackref,NtoStack,threshold,acTime,Offset)

% ProcessAc processes acoustic data, i.e., stacks the first 'NtoStackref' WFs to use as a
% reference signal and cross-correlate each WF with this reference to
% extract a change in time of flight. It also computes the amplitude of
% each waveform (peak-peak amplitude) and the root-mean-square amplitude of
% each WF.
% Waveforms can be stacked to improve SNR by choosing 'NtoStack' larger
% than 1.
% When choosing 'absoluteref', each stacked WF is cross-correlated with the
% first stacked WFs (made of the first NtoStackref waveforms).
% When choosing 'relativeref', each stacked WF is cross-correlated with the
% previous stacked WF. When choosing 'absoluteref', the window of interest
% is shifted by one sample when 'Timeshift' becomes larger than sampling
% time 'ts'.

% INPUTS
% Ac_path is the path where pXXXX.mat (containing acoustic settings) can be
% loaded
% run_ac_path is the path where the acoustic data can be loaded
% ts is the sampling time, typically it is the output 'ts_adjusted'
% returned by SyncAcData function
% totalnumberoffiles is the total number of acoustic files to be analyzed
% (output of SyncAcData)
% idxBeg is the idx corresponding to the beginning of the waveform to be
% analyzed
% idxEnd is the idx corresponding to the end of the waveform to be
% analyzed
% idx_TOF_0 is the arrival time of the reference waveform used to estimate
% a baseline time of flight (does not have to be precise for the routine to
% work properly)
% stepoptions is the ouput of the selectstep function
% reference is set either to 'absoluteref' or 'relativeref'. When using
% 'absoluteref', each WF is cross-correlated with the first NtoStackref WFs,
% where NtoStackref is the last input
% 'relativeref', each WF is cross-correlated with the previous WF.
% NtoStackref is the number of waveforms used to build a reference WF when
% using 'absoluteref'. When 'relativeref' is used, NtoStackref is set equal
% to NtoStack by default. 
% threshold is used to ignore noisy waveforms when using 'relativeref'. 
% TimeShift and RmsAmp values are kept only when MaxInter is above this
% threshold. threshold is comprised between -1 (no threshold) and 1 (max
% threshold). When 'absoluteref' is used, threshold is ignored.
% acTime is the acoustic time vector, typically the output from SyncAcData. 
% Offset is a number allowing one to offset waveforms
% corresponding to different channels

% OUTPUTS
% MaxInter is the maximum of intercorrelation vector
% MinInter is the minimum of intercorrelation vector
% TimeShift is the delay found between each WF with its reference WF
% (either absolute or relative reference). When the reference is relative,
% cumsum can be used to sum up the cumulative time differences.
% RmsAmp contains the rms amplitude of each WF
% TOF_0 is the reference time-of-flight that can be used to estimate a
% baseline velocity
% RmsAmpRef is the rms amplitude of the reference WF 
% fullWFref is the reference WF.
% acTime_new is the new acoustic time vector, which accounts for stacking.
% If NtoStack is 1, acTime is equal to the output acTime of SyncAcData. If
% NtoStack is larger than 1, it is modified to account for the stacking.
% integerShift indicates by how much the window is shifted for each
% waveform during the analysis (it will be useful only when the section at
% 303 is uncommented, otherwise it will be a matrix of zeros that you can
% ignore).

% if no offset is specified, set it to 0.
if nargin < 14
    Offset = 0;
end

if strcmp(reference,'relativeref')
    NtoStackref = NtoStack;
    warning('The number of reference waveforms (NtoStackref) when using ''relativeref'' is set equal to ''NtoStack''.')    
end

if strcmp(reference,'absoluteref') && (NtoStackref ~= 1)
    threshold = -1;
    warning('The threshold is set to -1 (i.e., no threshold) when using ''absoluteref''.')    
end

% totalnumberoffiles = 10; % uncomment to test the program with 10 acoustic files only

% acoustic parameters
acSettings = load(Ac_path);                                 % load acoustic settings
numSFpfile = acSettings.numFrames/2;                        % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;                           % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;                     % number of WF per file and per channel
numCHR = length(acSettings.channels2save);                   % number of receiving channels
numCHT = length(acSettings.channels2transmit);              % number of transmitting channels
WFlength = acSettings.Nsamples;                             % waveform length
fs = 1/ts;                                                  % acoustic sampling rate
clear acSettings

acN = totalnumberoffiles*numWFpfilepCH; % total number of WF per channel
acN = floor(acN/NtoStack/numCHT); % total number of stacked WF per channel and per transmitter

% modify acTime to account for stacking and the number of transmitters used
acPeriod = mean(diff(acTime));
acPeriod_new = NtoStack*acPeriod*numCHT;

acTime_new = 0:acPeriod_new:acPeriod_new*(acN-1);
acTime_newshifted = zeros(length(acTime_new),numCHT); % matrix containing one time vector per transmitter
for chnumt = 1:numCHT % shift by one acTime sample times number of WF to stack for each new transmitter
    acTime_newshifted(:,chnumt) = acTime_new + acTime(chnumt)*NtoStack; 
end

acTime_new = acTime_newshifted;

% time vector for each waveform
timeWF = (0:WFlength-1)*ts;

% apodization function
% pulse_length = idxEnd - idxBeg + 1;
% hanning_length = round(10*pulse_length/100);
% hann = [0; hanning(hanning_length-2); 0];
% window = ones(pulse_length,1);
% window(1:ceil(hanning_length/2)) = hann(1:ceil(hanning_length/2));
% window(pulse_length-floor(hanning_length/2)+1:pulse_length) = hann(ceil(hanning_length/2)+1:hanning_length);

fullWFref = zeros(WFlength,numCHR,numCHT);
RmsAmpRef = zeros(numCHR,numCHT);
AmpRef = zeros(numCHR,numCHT);

MaxInter = zeros(acN,numCHR,numCHT);
MinInter = zeros(acN,numCHR,numCHT);
tempTimeShift = zeros(acN,numCHR,numCHT); % temporary time shift
TimeShift = zeros(acN,numCHR,numCHT);
RmsAmp = zeros(acN,numCHR,numCHT);
Amp = zeros(acN,numCHR,numCHT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a reference waveform, based on the first NtoStackref WFs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if 'relativeref' is chosen, this reference WF will be used once only to
% be compared with the next one

kk = 0; % from 0 to NtoStackref*numCHT - 1
ii = 1; % number of files needed to build the reference WF
jj = 2; % from 1 to numWFpfilepCH (starts at 2 to avoid the very first WF)
chnumt = 1; % transmitter index
while kk < NtoStackref*numCHT % number of Reference WFs    
    if jj == 1 || (ii == 1 && jj == 2) % open new file if jj = 1 (except for the first file, open new file when jj = 2)        
        ACfilename = [run_ac_path num2str(ii) '.ac']; % only the first file is needed to extract the first 50 WF
        fid = fopen(ACfilename,'r');
        ACdata = fread(fid,'int16');
        fclose(fid);
        
        % reshape to get one column per channel
        ACdata = reshape(ACdata,[],numCHR,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
        ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
        ACdata = reshape(ACdata,[],numCHR,1); % WF vs Channel
    end        
    fullWFref(:,:,chnumt) = fullWFref(:,:,chnumt) + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % stack WFs
    if chnumt < numCHT % chnumt runs from 1 to numCHT
        chnumt = chnumt + 1;
    else
        chnumt = 1;
    end   
    if jj < numWFpfilepCH  % stay within the same file for the next run
        jj = jj + 1;
    else                    % use next file for the next run
        jj = 1;ii = ii + 1;
    end
    kk = kk + 1;
end

fullWFref = fullWFref/NtoStackref; % average
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
        plot(timeWF(idx_TOF_0),fullWFref(idx_TOF_0,chnumr,chnumt)-Offset*(idxoffset-1),'sk');
        legend('FullWF','WFtoAnalyze','Arrival');
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

TOF_0 = idx_TOF_0/fs; % reference time of flight
fprintf(['Reference Time of Flight is ' num2str(TOF_0) ' micro-seconds\n\n'])

if stepoptions.SHOW_WFref
    fprintf(['Stop routine to look at reference WF.\n' ...
        'Pick up idxBeg (beginning of window for analysis),\n' ...
        'idxEnd (end of window for analysis) and \n' ...
        'idx_TOF_0 (WF arrival).\n' ...
        'Once done, run the program again with step 5 or 6.\n'])
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute changes in time of flight, max of intercorrelation and RmsAmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 1; % file number
jj = 1; % WFnum within each file

WFlengthanalysis = idxEnd - idxBeg + 1;
WF = zeros(WFlengthanalysis,numCHR,numCHT);

idxBegRef = idxBeg;
idxEndRef = idxEnd;

idxBeg = idxBeg*ones(acN,numCHR,numCHT);
idxEnd = idxEnd*ones(acN,numCHR,numCHT);
WindowShift = zeros(acN,numCHR,numCHT); % window shifted to "follow" the waveform
integerShift = zeros(acN,numCHR,numCHT); % from one WF to the next one, it is shifted by integerShift samples
lastindexshifted = zeros(numCHR,numCHT); % last sample the window was shifted

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
    for kk = 1:NtoStack*numCHT      
        if jj == 1 % if first WF of the file
            fprintf(['File number ' num2str(ii) '.\n']) % display file number
            ACfilename = [run_ac_path num2str(ii) '.ac'];
            fid = fopen(ACfilename,'r');
            ACdata = fread(fid,'int16');
            fclose(fid);
            
            % reshape to get one column per channel
            ACdata = reshape(ACdata,[],numCHR,numSFpfile);   % 3D matrix with WF vs Channel vs number of SF
            ACdata = permute(ACdata,[1 3 2]);               % put Channel as the last dimension before reshaping
            ACdata = reshape(ACdata,[],numCHR,1);            % WF vs Channel
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
    
    % cross-correlate
    for chnumt = 1:numCHT
        for chnumr = 1:numCHR
            WF(:,chnumr,chnumt) = fullWF(idxBeg(hh,chnumr,chnumt):idxEnd(hh,chnumr,chnumt),chnumr,chnumt); % WF is only the part to be analyzed
            corr_signals = xcorr(WFref(:,chnumr,chnumt),WF(:,chnumr,chnumt),'coeff');
            [MaxInter(hh,chnumr,chnumt),tempTimeShift(hh,chnumr,chnumt),MinInter(hh,chnumr,chnumt),] = delay(corr_signals,ts);
            TimeShift(hh,chnumr,chnumt) = tempTimeShift(hh,chnumr,chnumt) - WindowShift(hh,chnumr,chnumt)*ts;
            RmsAmp(hh,chnumr,chnumt) = rms(WF(:,chnumr,chnumt));                        % RmsAmp of the waveform
            Amp(hh,chnumr,chnumt) = max(WF(:,chnumr,chnumt))-min(WF(:,chnumr,chnumt));  % Max Amp of the waveform            
        end
    end
    if hh/500==round(hh/500)    
        display(MaxInter(hh,:,:))
    end
    % display
    if (stepoptions.SHOW_allWF && ~stepoptions.SHOW_1WFperfile) || (~stepoptions.SHOW_allWF && stepoptions.SHOW_1WFperfile && jj == 1)
        idxoffset = 1;        
        for chnumt = 1:numCHT
            for chnumr = 1:numCHR                 
                h1(chnumr,chnumt) = plot(timeWF,fullWFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'r');hold on                
                h2(chnumr,chnumt) = plot(timeWF(idxBegRef:idxEndRef),WFref(:,chnumr,chnumt)-Offset*(idxoffset-1),'g');grid on;                
                h3(chnumr,chnumt) = plot(timeWF,fullWF(:,chnumr,chnumt)-Offset*(idxoffset-1),'b');
                h4(chnumr,chnumt) = plot(timeWF(idxBeg(hh,chnumr,chnumt):idxEnd(hh,chnumr,chnumt)),WF(:,chnumr,chnumt)-Offset*(idxoffset-1),'k');                                
                drawnow
                idxoffset = idxoffset + 1;                
            end
        end
        cla(gca); % clear former WFs without removing axis properties
    end   

%% when using absoluteref, you can uncomment below to shift the window and avoid large phase lag between WFs 
%     if strcmp(reference,'absoluteref')
%         if hh >= 10 && hh < acN
%             for chnumt = 1:numCHT
%                 for chnumr = 1:numCHR
%                     % shift window for the next loop if last 10 waveforms were larger than ts
%                     if sum(abs(tempTimeShift(hh-9:hh,chnumr,chnumt)) - ts > 0) == 10 && (hh - lastindexshifted(chnumr,chnumt) >= 10)
%                         lastindexshifted(chnumr,chnumt) = hh; % keep index for next time
%                         display(['T' num2str(chnumt) ' - R' num2str(chnumr)]);
%                         display(floor(mean(tempTimeShift(hh-9:hh,chnumr,chnumt))/ts))
%                         integerShift(hh,chnumr,chnumt) = floor(mean(tempTimeShift(hh-9:hh,chnumr,chnumt))/ts); % mean of the last 10 time shifts
%                         idxBeg(hh+1,chnumr,chnumt) = idxBeg(hh,chnumr,chnumt) - integerShift(hh,chnumr,chnumt);
%                         idxEnd(hh+1,chnumr,chnumt) = idxEnd(hh,chnumr,chnumt) - integerShift(hh,chnumr,chnumt);
%                         WindowShift(hh+1,chnumr,chnumt) = WindowShift(hh,chnumr,chnumt) - integerShift(hh,chnumr,chnumt);
%                         if idxBeg(hh+1,chnumr,chnumt) < 1
%                             idxBeg(hh+1,chnumr,chnumt) = 1;
%                             idxEnd(hh+1,chnumr,chnumt) = WFlengthanalysis;
%                             display('The window has shifted too fast towards early arrivals, probably because of some noisy waveforms.');
%                         end
%                         if idxEnd(hh+1,chnumr,chnumt) > WFlength
%                             idxBeg(hh+1,chnumr,chnumt) = WFlength-WFlengthanalysis+1;
%                             idxEnd(hh+1,chnumr,chnumt) = WFlength;
%                             display('The window has shifted too fast towards late arrivals, probably because of some noisy waveforms.');
%                         end
%                     else
%                         idxBeg(hh+1,chnumr,chnumt) = idxBeg(hh,chnumr,chnumt);
%                         idxEnd(hh+1,chnumr,chnumt) = idxEnd(hh,chnumr,chnumt);
%                         WindowShift(hh+1,chnumr,chnumt) = WindowShift(hh,chnumr,chnumt);
%                     end
%                 end
%             end
%         end
%     end    

    % choose reference WF for the next run when
    for chnumt = 1:numCHT
        for chnumr = 1:numCHR
            if strcmp(reference,'relativeref') && MaxInter(hh,chnumr,chnumt) >= threshold
                
                % current WF becomes reference WF for the next run if max
                % of intercorrelation is above threshold
                fullWFref(:,chnumr,chnumt) = fullWF(:,chnumr,chnumt);
                WFref(:,chnumr,chnumt) = WF(:,chnumr,chnumt);
                
            elseif strcmp(reference,'relativeref') && MaxInter(hh,chnumr,chnumt) < threshold
                
                % MaxInter is lower than threshold.
                % result of cross-correlation is replaced by NaN and we
                % don't keep the current WF as a reference for the next run
                TimeShift(hh,chnumr,chnumt) = NaN;
                RmsAmp(hh,chnumr,chnumt) = NaN;                
            end
        end
    end
end
% Amplitude found for first WF is always wrong. We assign it the amplitude
% of the second WF.
if NtoStack == 1
    RmsAmp(1,:,:) = RmsAmp(2,:,:);
    TimeShift(1,:,:) = TimeShift(2,:,:);
end

% take the opposite such that a positive TimeShift means later arrival.
TimeShift = -TimeShift;

% sum up all time differences when 'relativeref' is used
if strcmp(reference,'relativeref')
    idxNaN = isnan(TimeShift); % find NaNs
    TimeShift(idxNaN) = 0; % 0 instead of NaN to be able to use cumsum below
    
    % sum up all the delays and put back NaN into the final vector
    TimeShift = cumsum(TimeShift);
    TimeShift(idxNaN) = NaN;
end

end

