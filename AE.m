function [EvAmp,EvTime] = AE(Ac_path,run_ac_path,ts,acTime,totalnumberoffiles,ParamAE,showEvents)

% This function is used to analyze acoustic emission data

% Inputs
% Ac_path is the path where pXXXX.mat (containing acoustic settings) can be
% loaded
% run_ac_path is the path where the acoustic data can be loaded
% ts is the sampling time, typically it is the output 'ts_adjusted'
% returned by SyncAcData function
% acTime is the time vector synced with mechanical data. It is typically
% the output of SyncAcData function.

% Outputs
% The only output is a figure with displayed waveforms

% acoustic parameters
acSettings = load(Ac_path);                     % load acoustic settings
numSFpfile = acSettings.numFrames/2;            % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;               % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;         % number of WF per file and per channel
numCH = length(acSettings.channels2save);       % number of channels
WFlength = acSettings.Nsamples;                 % segment length
ts = ts/1e6;                                    % from microsec to sec
fs = 1/ts;                                      % acoustic sampling rate
clear acSettings

% time vector for each waveform
timeWF = (0:WFlength-1)'*ts;

AETime = NaN(WFlength*numWFpfilepCH,1);

Last10EventsPrevFileTime = zeros(10,numCH);
Last10EventsPrevFileAmp = zeros(10,numCH);

for ii = 1:totalnumberoffiles                    
        
    % build time vector for each file   
    for jj = 1:numWFpfilepCH
        AETime((jj-1)*WFlength+1:jj*WFlength) = acTime((ii-1)*numWFpfilepCH+jj) + timeWF;
    end    
    
    ACfilename = [run_ac_path num2str(ii) '.ac']; % only the first file is needed to extract the first 50 WF
    fid = fopen(ACfilename,'r');
    ACdata = fread(fid,'int16');
    fclose(fid);
    
    % reshape to get one column per channel
    ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
    ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
    ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel                   
    
    for kk = 1:numCH
    
        chname = ['ch' num2str(kk)];
        ACdata(:,kk) = ACdata(:,kk)-mean(ACdata(:,kk));                                    
        
        % extract envelope        
        Env = abs(hilbert(ACdata(:,kk)));
        % smooth envelope
        Env = smooth(Env,ParamAE.Nsmooth); % AETime, 100 is good
        [EvAmptmp,EvTimetmp] = findpeaks(Env,AETime,...
            'MinPeakDistance',ParamAE.minPkDist,'MinPeakHeight',ParamAE.minPkHeight); % ,'MinPeakProminence',20

        % display only 25ms of data per file
        if showEvents
            if ii == 1 && kk == 1, figure;   end
            offsetplot = 0; %1.2*max(max(ACdata(1:1:100000,:)));
            plot(AETime(1:1:100000),ACdata(1:1:100000,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on
            dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
            title(['Channel ' num2str(kk)]);
            plot(AETime(1:1:100000),Env(1:1:100000)+offsetplot*(kk-1));
            idx = find(EvTimetmp < AETime(100000),1,'last');
            plot(EvTimetmp(1:idx),EvAmptmp(1:idx)+offsetplot*(kk-1),'o');
        end
% %         % display full file (roughly 1s of data) when RD has already
% %         % been assigned
% %         if ii == 1 && kk == 1, figure;   end
% %         offsetplot = 0; %1.2*max(max(ACdata(1:10:end,:)));
% %         plot(AETime(1:10:end),ACdata(1:10:end,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on
% %         dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
% %         plot(AETime(1:10:end),Env(1:10:end)+offsetplot*(kk-1));
% %         title(['Channel ' num2str(kk)]);
% %         plot(EvTimetmp,EvAmptmp+offsetplot*(kk-1),'o');
                
        % The last 10 events of the previous file are added to the
        % beginning of the vector to analyze the first 10 events of the
        % current file
        if ii == 1 % The first 10 events of the first file will be kept
            EvAmptmp = [zeros(10,1); EvAmptmp];
            EvTimetmp = [zeros(10,1); EvTimetmp];            
        else 
            EvAmptmp = [Last10EventsPrevFileAmp(:,kk); EvAmptmp];
            EvTimetmp = [Last10EventsPrevFileTime(:,kk); EvTimetmp];
        end
        
        EvAmptmp2 = []; % events when 'ring-down' events are removed from EvAmptmp
        EvTimetmp2 = [];
        for jk = 11:length(EvTimetmp) % starts at 11 because the first ten events are now from the previous file (already counted)
            ThresholdPreviousEvents = NaN(10,1);
            for hg = 1:10 % threshold for the last 10 events
                ThresholdPreviousEvents(hg) = EvAmptmp(jk-hg)*exp(-(EvTimetmp(jk) - EvTimetmp(jk-hg))/ParamAE.threshold);
            end
            condition = EvAmptmp(jk) >= ThresholdPreviousEvents; % condition is a  vector of 10 elements having 0 or 1
            if sum(condition) == 10 % i.e. all 10 previous events have a lower ring-down threshold
                EvAmptmp2 = [EvAmptmp2; EvAmptmp(jk)];
                EvTimetmp2 = [EvTimetmp2; EvTimetmp(jk)];
            end            
        end
        Nevents = length(EvAmptmp2); % number of events
        fprintf(['We found ' num2str(Nevents) ' events for file ' num2str(ii) ', channel ' num2str(kk) '.\n']);        
               
        if showEvents
            idx = find(EvTimetmp2 < AETime(100000),1,'last');
            plot(EvTimetmp2(1:idx),EvAmptmp2(1:idx)+offsetplot*(kk-1),'*');hold off
            drawnow;
%             pause
        end
        
% %         plot(EvTimetmp2,EvAmptmp2+offsetplot*(kk-1),'*');hold off
% %         drawnow;
% %         pause
        
        if ii == 1 % first file
            EvAmp.(chname) = EvAmptmp2;
            EvTime.(chname) = EvTimetmp2;
        else
            EvAmp.(chname) = [EvAmp.(chname); EvAmptmp2];
            EvTime.(chname) = [EvTime.(chname); EvTimetmp2];            
        end        
    
        % keep last 10 events info from current file to analyze first 10
        % events of the next file
        if Nevents >= 10
            Last10EventsPrevFileTime(:,kk) = EvTimetmp2(end-9:end);
            Last10EventsPrevFileAmp(:,kk) = EvAmptmp2(end-9:end);
        else % less than 10 events during the current file            
            Last10EventsPrevFileTime(end-Nevents+1:end,kk) = EvTimetmp2;
            Last10EventsPrevFileAmp(end-Nevents+1:end,kk) = EvAmptmp2;
        end
    end   
    hold off    

end
%%
%     ACdatas = zeros(size(ACdata,1),size(ACdata,2)); % used for kurtosis


%         plot(AETime(1:1:100000),ACdatas(1:1:100000,kk)+offsetplot*(kk-1),'Color',[0.6 , 0.6 , 0.6]);hold on 
%         plot(AETime(1:1:100000),Kurtosis(1:1:100000)+offsetplot*(kk-1));

          % kurtosis
%         ACdatas(:,kk) = smooth(ACdata(:,kk),20); % 3 was ok... 
%         M = 100;
%         Kurtosis = zeros(WFlength*numWFpfilepCH,1);
%         for ll = M/2:WFlength*numWFpfilepCH-M/2                        
%             ACdatawin = ACdatas(ll-M/2+1:ll+M/2,kk);            
% %             Kurtosis(ll) = sum((ACdatawin - mean(ACdatawin)).^4)/((M-1)*std(ACdatawin,1)^4) - 3;            
%             Kurtosis(ll) = kurtosis(ACdatawin) - 3;
%         end  


end




