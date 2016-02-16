function [MaxInter,TimeShift,RmsAmp,Amp,TOF_0,RmsAmpRef,AmpRef,fullWFref] = ...
    ProcessAc(Ac_path,run_ac_path,ts,totalnumberoffiles,idxBeg,idxEnd,idx_TOF_0, ...
              stepoptions,reference,NtoStackref,NtoStack,threshold)

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
% a baseline time of flight
% stepoptions is the ouput of the selectstep function
% reference is set either to 'absoluteref' or 'relativeref'. When using
% 'absoluteref', each WF is cross-correlated with the first NtoStackref WFs,
% where NtoStackref is the last input
% 'relativeref', each WF is the previous WF.
% NtoStackref is the number of waveforms used to build a reference WF when
% using 'absoluteref'. When 'relativeref' is used, NtoStackref is set equal to NtoStack by
% default. 
% threshold is used to ignore noisy waveforms when using 'relativeref'.
% TimeShift and RmsAmp values are kept only when MaxInter is above this threshold.
% threshold is comprised between -1 (no threshold) and 1 (max threshold).
% When 'absoluteref' is used, threshold is ignored.

% OUTPUTS
% MaxInter is the maximum of intercorrelation vector
% TimeShift is the delay found between each WF with its reference WF
% (either absolute or relative reference). When the reference is relative,
% cumsum can be used to sum up the cumulative time differences.
% RmsAmp contains the rms amplitude of each WF
% TOF_0 is the reference time-of-flight that can be used to estimate a
% baseline velocity
% RmsAmpRef is the rms amplitude of the reference WF 
% fullWFref is the reference WF.

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
numCH = length(acSettings.channels2save);                   % number of channels
WFlength = acSettings.Nsamples;                             % waveform length
fs = 1/ts;                                         % acoustic sampling rate
clear acSettings

acN = totalnumberoffiles*numWFpfilepCH; % total number of WF per channel
acN = floor(acN/NtoStack); % total number of stacked WF per channel

% time vector for each waveform
timeWF = (0:WFlength-1)*ts;

fullWFref = zeros(WFlength,numCH);
RmsAmpRef = zeros(numCH);
AmpRef = zeros(numCH);

MaxInter = zeros(acN,numCH);
tempTimeShift = zeros(acN,numCH); % temporary time shift
TimeShift = zeros(acN,numCH);
RmsAmp = zeros(acN,numCH);
Amp = zeros(acN,numCH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a reference waveform, based on the first NtoStackref WFs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if 'relativeref' is chosen, this reference WF will be used once only to
% be compared with the next one

kk = 0; % from 0 to NtoStackref - 1
ii = 1; % number of files needed to build the reference WF
jj = 2; % from 1 to numWFpfilepCH (starts at 2 to avoid the very first WF)
while kk < NtoStackref % number of Reference WFs
    if jj == 1 || (ii == 1 && jj == 2) % open new file if jj = 1 (except for the first file, open new file when jj = 2)
        ACfilename = [run_ac_path num2str(ii) '.ac']; % only the first file is needed to extract the first 50 WF
        fid = fopen(ACfilename,'r');
        ACdata = fread(fid,'int16');
        fclose(fid);
        
        % reshape to get one column per channel
        ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
        ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
        ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel
    end    
    fullWFref = fullWFref + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % stack WFs
    if jj < numWFpfilepCH  % stay within the same file for the next run
        jj = jj + 1;
    else                    % use next file for the next run
        jj = 1;ii = ii + 1;
    end
    kk = kk + 1;
end

fullWFref = fullWFref/NtoStackref; % average
WFref = fullWFref(idxBeg:idxEnd,:); % part of the WF to be analyzed

for chnum = 1:numCH
    RmsAmpRef(chnum) = rms(WFref(:,chnum)); % RmsAmp of the reference waveform
    AmpRef(chnum) = max(WFref(:,chnum))-min(WFref(:,chnum)); % RmsAmp of the reference waveform
end

figure(3)
for chnum = 1:numCH
    subplot(numCH,1,chnum);plot(timeWF*1e6,fullWFref(:,chnum),'r');hold on;grid on
    subplot(numCH,1,chnum);plot(timeWF(idxBeg:idxEnd)*1e6,WFref(:,chnum),'g');
    subplot(numCH,1,chnum);plot(timeWF(idx_TOF_0)*1e6,fullWFref(idx_TOF_0,chnum),'sk');hold off
    legend('FullWF','WFtoAnalyze','Arrival');
    xlabel('Time ($\mu s$)','Interpreter','Latex');
    ylabel('Amplitude (bits)','Interpreter','Latex');

    dcmObj = datacursormode;
    set(dcmObj,'UpdateFcn',@GoodCursor);
    set(gca,'FontSize',16);
    drawnow
end

TOF_0 = idx_TOF_0/fs; % reference time of flight
fprintf(['Reference Time of Flight is ' num2str(TOF_0*1e6) ' micro-seconds\n\n'])

if stepoptions.SHOW_WFref
    fprintf(['Stop routine to look at reference WF.\n' ...
        'Pick up idxBeg (beginning of window for analysis),\n' ...
        'idxEnd (end of window for analysis) and \n' ...
        'idx_TOF_0 (WF arrival).\n' ...
        'Once done, run the program again with step 4 or 5.\n'])
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute changes in time of flight, max of intercorrelation and RmsAmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fullWF = zeros(WFlength,numCH);

ii = 1; % file number
jj = 1; % WFnum within each file

WindowShift = 0;

idxBegRef = idxBeg;
idxEndRef = idxEnd;

for hh = 1:acN % from 1 to the total number of stacked waveforms
    % stack WFs
    for kk = 1:NtoStack        
        if jj == 1 % if first WF of the file
            fprintf(['File number ' num2str(ii) '.\n']) % display file number
            ACfilename = [run_ac_path num2str(ii) '.ac'];
            fid = fopen(ACfilename,'r');
            ACdata = fread(fid,'int16');
            fclose(fid);
            
            % reshape to get one column per channel
            ACdata = reshape(ACdata,[],numCH,numSFpfile);   % 3D matrix with WF vs Channel vs number of SF
            ACdata = permute(ACdata,[1 3 2]);               % put Channel as the last dimension before reshaping
            ACdata = reshape(ACdata,[],numCH,1);            % WF vs Channel
        end
        fullWF = fullWF + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % read data
        
        if jj < numWFpfilepCH   % stay within the same file for the next run
            jj = jj + 1;
        else                    % use next file for the next run
            jj = 1;ii = ii + 1;
        end
    end
    fullWF = fullWF/NtoStack; % stacked WF        
    WF = fullWF(idxBeg:idxEnd,:); % WF is only the part to be analyzed
     
    % cross-correlate
    for chnum = 1:numCH
        corr_signals = xcorr(WFref(:,chnum),WF(:,chnum),'coeff');
        [MaxInter(hh,chnum),tempTimeShift(hh,chnum)] = delay(corr_signals,ts);
        TimeShift(hh,chnum) = tempTimeShift(hh,chnum) - WindowShift*ts;
        RmsAmp(hh,chnum) = rms(WF(:,chnum)); % RmsAmp of the waveform
        Amp(hh,chnum) = max(WF(:,chnum))-min(WF(:,chnum)); % Max Amp of the waveform
    end
    % display
    if (stepoptions.SHOW_allWF && ~stepoptions.SHOW_1WFperfile) || (~stepoptions.SHOW_allWF && stepoptions.SHOW_1WFperfile && jj == 1)
        figure(3)
        for chnum = 1:numCH
            subplot(numCH,1,chnum);plot(timeWF*1e6,fullWFref(:,chnum),'r');hold on
            subplot(numCH,1,chnum);plot(timeWF(idxBegRef:idxEndRef)*1e6,WFref(:,chnum),'g');grid on;
        end
        figure(3)
        for chnum = 1:numCH
            subplot(numCH,1,chnum);plot(timeWF*1e6,fullWF(:,chnum),'b');
            subplot(numCH,1,chnum);plot(timeWF(idxBeg:idxEnd)*1e6,WF(:,chnum),'k');hold off
            drawnow
        end
    end    
    % shift window for the next loop
    if strcmp(reference,'absoluteref') && abs(tempTimeShift(hh,1)) > ts %&& MaxInter(hh,1) > 0.5
        integerShift = floor(tempTimeShift(hh,1)/ts); % first channel
        idxBeg = idxBeg - integerShift;
        idxEnd = idxEnd - integerShift;
        WindowShift = WindowShift - integerShift;       
        if idxBeg < 1 
            error('The window has shifted too fast towards early arrivals, probably because of some noisy waveforms. Consider increasing the number of waveforms to stack to reduce noise');        
        end
    end    
    % choose reference WF for the next run
    for chnum = 1:numCH
        if strcmp(reference,'relativeref') && MaxInter(hh,chnum) >= threshold
            % current WF becomes reference WF for the next run if max
            % of intercorrelation is above threshold
            fullWFref(:,chnum) = fullWF(:,chnum);
            WFref(:,chnum) = WF(:,chnum);
        elseif strcmp(reference,'relativeref') && MaxInter(hh,chnum) < threshold
            % MaxInter is lower than threshold.
            % result of cross-correlation is replaced by NaN and we
            % don't keep the current WF as a reference for the next run
            TimeShift(hh,chnum) = NaN;
            RmsAmp(hh,chnum) = NaN;
        end
    end           
end

% Amplitude found for first WF is always wrong. We assign it the amplitude
% of the second WF.
if NtoStack == 1
    RmsAmp(1,:) = RmsAmp(2,:);
    TimeShift(1,:) = TimeShift(2,:);
end

% take the opposite such that a positive TimeShift means later arrival.
TimeShift = -TimeShift;

end





