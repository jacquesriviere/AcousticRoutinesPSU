function [MaxInter,TimeShift,RmsAmp,TOF_0,RmsAmpRef,fullWFref] = ProcessAc(Ac_path,run_ac_path,ts_adjusted,totalnumberoffiles,idxBeg,idxEnd,idx_TOF_0,stepoptions)

% ProcessAc processes acoustic data, i.e., use the first 50 WFs as a
% reference signal and cross-correlate each WF with this reference to
% extract a change in time of flight. It also compute the root-mean-square
% amplitude of each WF.

% acoustic parameters
acSettings = load(Ac_path);                                 % load acoustic settings
numSFpfile = acSettings.numFrames/2;                        % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;                           % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;                     % number of WF per file and per channel
numCH = length(acSettings.channels2save);                   % number of channels
WFlength = acSettings.Nsamples;                             % waveform length
fs = 1/ts_adjusted;                                         % acoustic sampling rate
clear acSettings

acN = totalnumberoffiles*numWFpfilepCH; % total number of WF per channel

fullWFref = zeros(WFlength,numCH);
RmsAmpRef = zeros(numCH);
MaxInter = zeros(acN,numCH);
TimeShift = zeros(acN,numCH);
RmsAmp = zeros(acN,numCH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build a reference waveform, based on the first 50 WF.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
ACfilename = [run_ac_path num2str(1) '.ac']; % only the first file is needed to extract the first 50 WF
fid = fopen(ACfilename,'r');
ACdata = fread(fid,'int16');
fclose(fid);

% reshape to get one column per channel
ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel

for jj = 2:50 % do not use the very first waveform because the amplitude voltage is not correct
    fullWFref = fullWFref + ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % read data
end
fullWFref = fullWFref/49;
WFref = fullWFref(idxBeg:idxEnd,:);

for chnum = 1:numCH
    RmsAmpRef(chnum) = rms(WFref(:,chnum)); % RmsAmp of the waveform
end

figure(3)
for chnum = 1:numCH
    subplot(numCH,1,chnum);plot(fullWFref(:,chnum),'r');hold on;grid on
    subplot(numCH,1,chnum);plot(idxBeg:idxEnd,WFref(:,chnum),'g');
    subplot(numCH,1,chnum);plot(idx_TOF_0,fullWFref(idx_TOF_0,chnum),'sk');hold off
    legend('FullWF','WFtoAnalyze','Arrival');
    drawnow
end

if stepoptions.SHOW_WFref
    fprintf(['Stop routine to look at reference WF.\n' ...
        'Pick up idxBeg (beginning of window for analysis),\n' ...
        'idxEnd (end of window for analysis) and \n' ...
        'idx_TOF_0 (WF arrival).\n' ...
        'Once done, run the program again with step 4 or 5.\n'])
    %         error('stop routine to look at reference WF') % stop routine to pick up idx_TOF_0
    return
end

TOF_0 = idx_TOF_0/fs; % reference time of flight
fprintf(['Reference Time of Flight is ' num2str(TOF_0*1e6) ' micro-seconds\n\n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute changes in time of flight, max of intercorrelation and RmsAmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hh = 1;
for ii = 1:totalnumberoffiles
    fprintf(['File number ' num2str(ii) '.\n'])
    ACfilename = [run_ac_path num2str(ii) '.ac'];
    fid = fopen(ACfilename,'r');
    ACdata = fread(fid,'int16');
    fclose(fid);
    
    % reshape to get one column per channel    
    ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF   
    ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
    ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel
            
    for jj = 1:numWFpfilepCH        
        fullWF = ACdata(WFlength*(jj-1)+1:WFlength*jj,:); % read data
        WF = fullWF(idxBeg:idxEnd,:); % WF is only the part to be analyzed
      
        for chnum = 1:numCH
            corr_signals = xcorr(WFref(:,chnum),WF(:,chnum),'coeff');
            [MaxInter(hh,chnum),TimeShift(hh,chnum)] = delay(corr_signals,ts_adjusted);            
            RmsAmp(hh,chnum) = rms(WF(:,chnum)); % RmsAmp of the waveform            
        end
        hh = hh + 1;        
        if (stepoptions.SHOW_allWF && ~stepoptions.SHOW_1WFperfile) || (~stepoptions.SHOW_allWF && stepoptions.SHOW_1WFperfile && jj == 1)            
            figure(3)
            for chnum = 1:numCH
                subplot(numCH,1,chnum);plot(fullWFref(:,chnum),'r');hold on
                subplot(numCH,1,chnum);plot(idxBeg:idxEnd,WFref(:,chnum),'g');grid on;
            end                        
            figure(3)
            for chnum = 1:numCH
                subplot(numCH,1,chnum);plot(fullWF(:,chnum),'b');
                subplot(numCH,1,chnum);plot(idxBeg:idxEnd,WF(:,chnum),'k');hold off
                drawnow
            end            
        end
    end    
end

% Amplitude found for first WF is always wrong. We assign it the amplitude
% of the second WF.
RmsAmp(1) = RmsAmp(2);

% take the opposite such that a positive TimeShift means later arrival.
TimeShift = -TimeShift;


end

