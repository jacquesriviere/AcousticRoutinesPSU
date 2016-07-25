function RD = setupAEdetection(Ac_path,run_ac_path,ts,acTime,AtWhichTimes,ParamAE)

% This function is used to setup parameters in order to properly detect AE events
% and time arrivels when using the main computeAE function

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

N = length(AtWhichTimes);
leg = zeros(N,1);
rangeTimes = zeros(N,2);
legendmatrix = cell(N,1);
color = lines(N);
hh = 1; % index for legend
indexlegend = NaN(N,1);


% leg = zeros(totalnumberoffiles,1);
% 
% legendmatrix = cell(totalnumberoffiles,1);
% color = lines(totalnumberoffiles);
% hh = 1; % index for legend
% indexlegend = NaN(totalnumberoffiles,1);
AETime = NaN(WFlength*numWFpfilepCH,1);
rate = NaN(3*numCH,length(1000:100:1100)); % 3 for 3 events per figure

for ii = 1:N

    idxAcTime = find(acTime > AtWhichTimes(ii),1,'first');    
    if acTime(1) > AtWhichTimes(ii)      
        error(['The first value of ''AtWhichTimes'' is too small. No acoustic data at time ' num2str(AtWhichTimes(1)) 's.']);
    elseif acTime(end) < AtWhichTimes(ii)
        error(['The last value of ''AtWhichTimes'' is too large. No acoustic data at time ' num2str(AtWhichTimes(end)) 's.']);
    end             
    
    filenumber = ceil(idxAcTime/numWFpfilepCH); % file number for the first WF to stack    
        
    % build time vector for each file   
    for jj = 1:numWFpfilepCH
        AETime((jj-1)*WFlength+1:jj*WFlength) = acTime((filenumber-1)*numWFpfilepCH+jj) + timeWF;
    end    
    
    ACfilename = [run_ac_path num2str(filenumber) '.ac'];
    fid = fopen(ACfilename,'r');
    ACdata = fread(fid,'int16');   
    fclose(fid);
    
    % reshape to get one column per channel
    ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
    ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
    ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel                   
    
%   ACdatas = zeros(size(ACdata,1),size(ACdata,2)); % used for kurtosis
    for kk = 1:numCH
        
%       chname = ['ch' num2str(kk)];
        ACdata(:,kk) = ACdata(:,kk)-mean(ACdata(:,kk)); % remove mean values             
        
        % extract envelope        
        Env = abs(hilbert(ACdata(:,kk)));
        % smooth envelope
        Env = smooth(Env,ParamAE.Nsmooth); % AETime, 100 is good
        [EvAmptmp,EvTimetmp] = findpeaks(Env,AETime,...
            'MinPeakDistance',ParamAE.minPkDist,'MinPeakHeight',ParamAE.minPkHeight); % ,'MinPeakProminence',20                                     
        
        if isnan(ParamAE.threshold) % pick 6 coordinates
            % plot all only ~25ms of data per file when picking ring-down times            
            if ii == 1 && kk == 1, figure;   end
            idx = find(EvTimetmp < AETime(100000),1,'last');
            offsetplot = 0; %1.2*max(max(ACdata(1:1:100000,:)));
            plot(AETime(1:1:100000),ACdata(1:1:100000,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on
            plot(AETime(1:1:100000),Env(1:1:100000)+offsetplot*(kk-1));           
            title(['Channel ' num2str(kk)]);
            dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
            plot(EvTimetmp(1:idx),EvAmptmp(1:idx)+offsetplot*(kk-1),'o');                        
            drawnow; % pause
            
            [t,A] = ginput(6); % pick 6 coordinates, i.e. 3 events and their ring-downs
            t0 = t(1:2:end);
            A0 = A(1:2:end);
            t =  t(2:2:end);
            A = A(2:2:end);
            
            ratetmp = -log(A./A0)./(t - t0);
            [A0_sorted,II] = sort(A0);
            ratetmp = ratetmp(II); % sort rate based on increasing amplitude
            rate(3*(kk-1)+1:3*kk,ii) = ratetmp;
            
        else % if threshold value has been assigned, display events which have been kept
            
            EvAmptmp2 = EvAmptmp(1:10);
            EvTimetmp2 = EvTimetmp(1:10);
            for jk = 11:length(EvTimetmp)
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
            % display full file (roughly 1s of data) when RD has already
            % been assigned
            if ii == 1 && kk == 1, figure;   end
            offsetplot = 0; %1.2*max(max(ACdata(1:10:end,:)));
            plot(AETime(1:10:end),ACdata(1:10:end,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on
            plot(AETime(1:10:end),Env(1:10:end)+offsetplot*(kk-1));            
            title(['Channel ' num2str(kk)]);
            dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
            plot(EvTimetmp,EvAmptmp+offsetplot*(kk-1),'o');
            plot(EvTimetmp2,EvAmptmp2+offsetplot*(kk-1),'*');hold off
            drawnow;
            pause
%             uncomment to display the first 100000 points of the file
%             if ii == 1 && kk == 1, figure;   end
%             offsetplot = 0; %1.2*max(max(ACdata(1:10:end,:)));
%             idx = find(EvTimetmp < AETime(100000),1,'last');
%             idx2 = find(EvTimetmp2 < AETime(100000),1,'last');
%             plot(AETime(1:100000),ACdata(1:100000,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on
%             plot(AETime(1:100000),Env(1:100000)+offsetplot*(kk-1));            
%             title(['Channel ' num2str(kk)]);
%             dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
%             plot(EvTimetmp(1:idx),EvAmptmp(1:idx)+offsetplot*(kk-1),'o');
%             plot(EvTimetmp2(1:idx2),EvAmptmp2(1:idx2)+offsetplot*(kk-1),'*');hold off
%             drawnow;
%             pause
        
        
        end          
    end

    hold off    

end

if isnan(ParamAE.threshold)
    RDtime = 1./rate; % ring-down time in s   
    figure;
    plot(RDtime*1e3);ylabel('Ring-down time (ms)');xlabel('Event Number');
    
    RDtime = reshape(RDtime,1,[]);
    meanRDtime = mean(RDtime);
    stdRDtime = std(RDtime,1);
    threshold = meanRDtime + 4*stdRDtime;
    RD.meanRDtime = meanRDtime;
    RD.stdRDtime = stdRDtime;
    RD.threshold = threshold;
    fprintf(['Threshold for ring-down is found to be ' num2str(RD.threshold*1e3) ' ms.\n' ...
        'Assign this value (in seconds) to ''ParamAE.threshold'' and run the program again']);
else
    RD.meanRDtime = NaN;
    RD.stdRDtime = NaN;
    RD.threshold = ParamAE.threshold;
    fprintf(['Threshold for ring-down was formerly found to be ' num2str(RD.threshold*1e3) ' ms.\n']);
end
% display legend for one channel only
% legend(leg(indexlegend),legendmatrix)
end











%%
          % kurtosis
%         ACdatas(:,kk) = smooth(ACdata(:,kk),20); % 3 was ok... 
%         M = 100;
%         Kurtosis = zeros(WFlength*numWFpfilepCH,1);
%         for ll = M/2:WFlength*numWFpfilepCH-M/2                        
%             ACdatawin = ACdatas(ll-M/2+1:ll+M/2,kk);            
% %             Kurtosis(ll) = sum((ACdatawin - mean(ACdatawin)).^4)/((M-1)*std(ACdatawin,1)^4) - 3;            
%             Kurtosis(ll) = kurtosis(ACdatawin) - 3;
%         end
        
        
%         if ii == 1 % first file
%             EvAmp.(chname) = EvAmptmp;
%             EvTime.(chname) = EvTimetmp;
%         else
%             EvAmp.(chname) = [EvAmp.(chname); EvAmptmp];
%             EvTime.(chname) = [EvTime.(chname); EvTimetmp];            
%         end   

        % display full file (roughly 1s of data)
%         figure(4);
%         offsetplot = 1.2*max(max(ACdata(1:10:end,:)));
%         plot(AETime(1:10:end),ACdata(1:10:end,kk)+offsetplot*(kk-1),'Color',[0.2 , 0.2 , 0.2]);hold on  
%         plot(AETime(1:10:end),Env(1:10:end)+offsetplot*(kk-1));
%         dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
%         plot(EvTimetmp,EvAmptmp+offsetplot*(kk-1),'o');
%         drawnow; pause
        
        % display only 25ms of data per file
%         if ii == 1, figure;   end         
%         idx = find(EvTimetmp < AETime(100000),1,'last');
%         offsetplot = 0; %1.2*max(max(ACdata(1:1:100000,:)));        
%         plot(AETime(1:1:100000),ACdata(1:1:100000,kk)+offsetplot*(kk-1),'Color',[0.8 , 0.8 , 0.8]);hold on 
% %         plot(AETime(1:1:100000),ACdatas(1:1:100000,kk)+offsetplot*(kk-1),'Color',[0.6 , 0.6 , 0.6]);hold on 
% %         plot(AETime(1:1:100000),Kurtosis(1:1:100000)+offsetplot*(kk-1));
%         plot(AETime(1:1:100000),Env(1:1:100000)+offsetplot*(kk-1));
%         dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
%         plot(EvTimetmp(1:idx),EvAmptmp(1:idx)+offsetplot*(kk-1),'o');hold off
%         drawnow; % pause  
