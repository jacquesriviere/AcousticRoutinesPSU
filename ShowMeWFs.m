function ShowMeWFs(WF_path,AcSettingsfile,SyncFile,AtWhichTimes,NtoStack,Offset,WhichTransTomoMode)

% This function is used to display waveforms at particular times during the
% run, with the options of stacking them. It is typically used to pick
% arrival times by hand, or decide how many waveforms need to be stacked
% before further processing

% Inputs

% WF_path is the path where the acoustic data can be loaded

% AcSettingsfile is the name of the file containing the acoustic settings
% chosen when recording on Verasonics

% SyncFile is the result of the pXXXX_sa routine, providing the acTime vector
% and the adjusted sampling time ts returned by SyncAcData function

% AtWhichTimes is a scalar or vector of time values at which you wish to
% display the waveform

% NtoStack is the number of waveforms to stack, centered around each value
% of AtWhichTimes. The time range is indicated on the figure.

% Offset is the number allowing one to offset waveforms corresponding to
% different channels (put 0 if you don't want to offset waveforms)

% WhichTransTomoMode: When using tomography mode (i.e. transmit first
% emitter, record all receivers; transmit next emittor, record all
% receivers), ShowMeWFs will display waveforms recorded when transmitter
% "WhichTransTomoMode" was used. For instance, if transmitters 1, 2, 3 were
% P, SV and SH, WhichTransTomoMode = 1 will display WFs from all receivers
% when P-wave transmitter was active. If plane wave mode is used (i.e.
% transmit all emittors, record all receivers), do not assign any input and
% WhichTransTomoMode will be set to 1.

% Outputs
% The only output is a figure with displayed waveforms

if nargin < 7
    WhichTransTomoMode = 1;
end

% acoustic parameters
acSettings = load(AcSettingsfile);          % load acoustic settings
numSFpfile = acSettings.numFrames/2;        % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;           % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;     % number of WF per file and per channel
numCHR = length(acSettings.channels2save);   % number of channels
numCHT = length(acSettings.channels2transmit);   % number of channels

load(SyncFile); % load sync data

if WhichTransTomoMode > numCHT    
    error(['You chose to display transmitter #' num2str(WhichTransTomoMode) '. Please choose a transmitter between 1 and ' num2str(numCHT) '.']);
end

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
color = lines(N);
hh = 1; % index for legend
indexlegend = NaN(N,1);
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
    
    filenumber = ceil(idxAcTimeVec(1)/numWFpfilepCH); % file number for the first WF to stack
    idxWFwithinfile = mod(idxAcTimeVec(1),numWFpfilepCH); % closest WF within that file (but not necessarily corresponding to the good transmitter)
    
    % find the waveform which corresponds to transmitter WhichTransTomoMode
    ShiftTrans = mod(idxWFwithinfile,numCHT);
    if ShiftTrans == 0
        ShiftTrans = numCHT;  % i.e. no shift 
    end
    idxWFwithinfile = idxWFwithinfile + (WhichTransTomoMode - ShiftTrans); % closest WF within that file corresponding to the good trnasmitter         
    fullWFref = zeros(WFlength,numCHR);
    
    for kk = 1:NtoStack % number of WFs to stack                                
        
        if idxWFwithinfile <= numCHT || kk == 1 % open new file if idxWFwithinfile is 1 or if it's the first WF to stack
            ACdata = LoadAcFile(WF_path,filenumber,numCHR,numSFpfile);
        end
        
        fullWFref = fullWFref + ACdata(WFlength*(idxWFwithinfile-1)+1:WFlength*idxWFwithinfile,:); % stack WFs
        if idxWFwithinfile <= numWFpfilepCH - numCHT % stay within the same file for the next loop
            idxWFwithinfile = idxWFwithinfile + numCHT; % update to the next WF corresponding to the same transmitter 
        else                                % use next file for the next loop
            filenumber = filenumber + 1; % go to next file
            idxWFwithinfile = WhichTransTomoMode; % start in the next file at WFs corresponding to transmitter "WhichTransTomoMode"
        end        
    end   
    fullWFref = fullWFref/NtoStack;    
    if ii == 1 % create new figure the first time only        
        figure;
    end
    indexlegend(ii) = hh;
    for chnum = 1:numCHR
%         subplot(numCHR,1,chnum);
        leg(hh) = plot(timeWF,fullWFref(:,chnum)-Offset*(chnum-1),'linewidth',2,'col',color(ii,:));hold on;grid on
%         ylabel('Amplitude (bits)','Interpreter','Latex');
        set(gca,'FontSize',16);
        set(gca,'YTickLabel',[])
        dcmObj = datacursormode;
        set(dcmObj,'UpdateFcn',@GoodCursor);
        drawnow        
        hh = hh + 1; % index for legend        
    end
    xlabel('Time ($\mu s$)','Interpreter','Latex'); 
    title(['Transmitter ' num2str(WhichTransTomoMode)])
    legendmatrix{ii}=['Stacked over time range ' num2str(rangeTimes(ii,1),'%.2f') '-' num2str(rangeTimes(ii,2),'%.2f') ' s.'];
end
% display legend for one channel only
legend(leg(indexlegend),legendmatrix)
end




