function ShowMeWFs(Ac_path,run_ac_path,ts,acTime,AtWhichTimes,NtoStack)

% This function is used to display waveforms at particular times during the
% run, with the options of stacking them. It is typically used to pick
% arrival times by hand, or decide how many waveforms need to be stacked
% before further processing (cross-correlation and so on)

% Inputs
% Ac_path is the path where pXXXX.mat (containing acoustic settings) can be
% loaded
% run_ac_path is the path where the acoustic data can be loaded
% ts is the sampling time, typically it is the output 'ts_adjusted'
% returned by SyncAcData function
% acTime is the time vector synced with mechanical data. It is typically
% the output of SyncAcData function.
% AtWhichTimes is a scalar or vector of time values at which you wish to
% display the waveform
% NtoStack is the number of waveforms to stack, centered around each value
% of AtWhichTimes. The time range is indicated on the figure.

% Outputs
% The only output is a figure with displayed waveforms



% acoustic parameters
acSettings = load(Ac_path);                 % load acoustic settings
numSFpfile = acSettings.numFrames/2;        % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;           % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;     % number of WF per file and per channel
numCH = length(acSettings.channels2save);   % number of channels
WFlength = acSettings.Nsamples;             % waveform length
fs = 1/ts;                                  % acoustic sampling rate
clear acSettings

% time vector for each waveform
timeWF = (0:WFlength-1)*ts;

% number of WFs to show
N = length(AtWhichTimes);
leg = zeros(N,1);
rangeTimes = zeros(N,2);
legendmatrix = cell(N,1);
for ii = 1:N
    idxAcTime = find(acTime > AtWhichTimes(ii),1,'first');      
    idxAcTimeVec = idxAcTime-ceil(NtoStack/2)+1:idxAcTime+floor(NtoStack/2); % vector of indexes centered around 'idxAcTime'    
    if idxAcTimeVec(1) < 1
        error(['Either the first value of ''AtWhichTimes'' is too small or ''NtoStack'' is too large. Can not stack ' ...
            num2str(NtoStack) ' waveforms centered around ' num2str(AtWhichTimes(1)) ' s.']);
    elseif idxAcTimeVec(end) > length(acTime)
        error(['Either the last value of ''AtWhichTimes'' is too large or ''NtoStack'' is too large. Can not stack ' ...
            num2str(NtoStack) ' waveforms centered around ' num2str(AtWhichTimes(end)) ' s.']);
    end
    
    rangeTimes(ii,:) = [acTime(idxAcTimeVec(1)) acTime(idxAcTimeVec(end))];     
    
    filenumber = floor(idxAcTimeVec(1)/numWFpfilepCH); % file number for the first WF to stack

    idxWFwithinfile = mod(idxAcTime,numWFpfilepCH); 
    fullWFref = zeros(WFlength,numCH);
    for kk = 1:NtoStack % number of WFs to stack                                
        
        if idxWFwithinfile == 1 || kk == 1 % open new file if idxWFwithinfile is 1 or if it's the first WF to stack
            ACfilename = [run_ac_path num2str(filenumber) '.ac']; % only the first file is needed to extract the first 50 WF
            fid = fopen(ACfilename,'r');
            ACdata = fread(fid,'int16');
            fclose(fid);
            
            % reshape to get one column per channel
            ACdata = reshape(ACdata,[],numCH,numSFpfile); % 3D matrix with WF vs Channel vs number of SF
            ACdata = permute(ACdata,[1 3 2]); % put Channel as the last dimension before reshaping
            ACdata = reshape(ACdata,[],numCH,1); % WF vs Channel
        end
        
        fullWFref = fullWFref + ACdata(WFlength*(idxWFwithinfile-1)+1:WFlength*idxWFwithinfile,:); % stack WFs
        if idxWFwithinfile < numWFpfilepCH  % stay within the same file for the next loop
            idxWFwithinfile = idxWFwithinfile + 1;
        else                                % use next file for the next loop
            idxWFwithinfile = 1;filenumber = filenumber + 1;
        end        
    end   
    fullWFref = fullWFref/NtoStack;
    if ii == 1 % create new figure the first time only
        figure;
    end
    for chnum = 1:numCH
        subplot(numCH,1,chnum);
        leg(ii) = plot(timeWF*1e6,fullWFref(:,chnum));hold on;grid on
        ylabel('Amplitude (bits)','Interpreter','Latex');
        set(gca,'FontSize',16);
        dcmObj = datacursormode;
        set(dcmObj,'UpdateFcn',@GoodCursor);
        drawnow
    end
    xlabel('Time ($\mu s$)','Interpreter','Latex'); 
    legendmatrix{ii}=['Stacked over time range ' num2str(rangeTimes(ii,1),'%.2f') '-' num2str(rangeTimes(ii,2),'%.2f') ' s.'];
end
legend(legendmatrix)

end




