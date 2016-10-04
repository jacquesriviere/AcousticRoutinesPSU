function [acTime,acPeriod_adjusted,ts_adjusted,TotalNumberOfFiles] = SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2)
% SyncAcData uses acoustic settings (mat file from Verasonics) and indexes
% provided by user to sync mechanical and acoustical data. It accounts for
% small variations in sampling frequency/period (of the order of
% picoseconds). It uses the biax recorder (mechanical data) as a
% reference and assume that the sampling frequency for acoustical data (and
% therefore acoustic pulsing rate) is slightly off.

% The main output is a time vector for acoustic data. It also outputs an
% adjusted acoustic sampling rate and an adjusted acoustic pulsing rate.

% INPUTS
% AcSettingsfile is the path where pXXXX.mat (containing acoustic settings) can be
% loaded
% Time is the time vector from biax data
% Sync is the time vector from biax data
% idxft is the index corresponding to the first trigger of the run
% idxlt is the index corresponding to the last trigger of the run
% idxref1 and idxref2 are two large triggers chosen towards the beginning
% and the end of the run

% OUTPUTS
% acTime is the time vector for acoustic data obtained after synchronization
% acPeriod_adjusted is the time between two consecutive acoustic pulses in microsec (found after adjusting using biax recorder as a reference)
% ts_adjusted is the sampling period in microsec (i.e. 1/fs) found after adjusting using biax recorder as a reference
% TotalNumberOfFiles is the total number of acoustic files to be analyzed 

% acoustic parameters
acSettings = load(AcSettingsfile);                     % load acoustic settings
acPeriod = acSettings.SeqControl(1).argument;   % time btw pulses in microsec (SeqControl(1).argument = timeBetweenpulses)
acRate = 1e6/(acPeriod);                        % in Hz (number of WF per second)
numSFpfile = acSettings.numFrames/2;            % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;               % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;         % number of WF per file and per channel
fs = acSettings.samplingFreq;                   % acoustic sampling rate in MHz
clear acSettings

ts = 1/fs; % sampling time in microsec

MAX = max(Sync(idxft:idxlt));
MIN = min(Sync(idxft:idxlt));
AMP = MAX-MIN;
Sync(idxft:idxlt) = (Sync(idxft:idxlt) - mean(Sync(idxft:idxlt)))/AMP;
Sync(idxft:idxlt) = Sync(idxft:idxlt) - max(Sync(idxft:idxlt));

figure;
plot(Time(idxft:idxlt),Sync(idxft:idxlt));hold on
xlabel('Time (s)');ylabel('Normalized Sync');
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);

% theoretical time between two consecutive triggers
t_btw_trig = numWFpfilepCH*acPeriod/1e6; % time in sec

% build a trigger time vector using the large trigger chosen above as a reference

% number of triggers between the two reference triggers:
% this number should be an integer if no mismatch occurs between mechanical and acoustical recorders.
% In the following, we take the biax recorder as the reference and adjust
% the acoustical data to match the mechanical data.
Nspans = (Time(idxref2) - Time(idxref1))/t_btw_trig;

% actual time between two consecutive triggers (using biax recorder as a reference)
actual_t_btw_trig = (Time(idxref2) - Time(idxref1))/round(Nspans); % sec   (round instead of floor here...)
format long
mismatchPRF = (actual_t_btw_trig - t_btw_trig)*1e6/numWFpfilepCH;
acPeriod_adjusted = acPeriod+mismatchPRF; % microsec
acRate_adjusted = 1e6/(acPeriod_adjusted); % Hz

fprintf(['In theory, acoustic pulses are sent every ' num2str(acPeriod/1e3,15) ' ms.\n' ...
         'Using biax recorder as a reference, we find ' num2str(acPeriod_adjusted/1e3,15) ' ms.\n' ...
         'Mismatch for the acoustic pulsing rate is ' num2str(mismatchPRF*1e3) ' ns.\n\n']);

mismatchts = (acRate/1e6)/fs*mismatchPRF; % microsec * MHz/MHz = microsec
ts_adjusted = ts+mismatchts; % adjusted sampling period (microsec)
fs_adjusted = 1/ts_adjusted; % adjusted sampling frequency (MHz)

fprintf(['Said differently, sampling frequency is ' num2str(fs) ' MHz in theory.\n' ...
         'Using biax recorder as a reference, it is adjusted to ' num2str(fs_adjusted,10) ' MHz.\n' ...
         'Mismatch for sampling time is ' num2str(mismatchts*1e6) ' ps.\n\n']);
     
% idxref2 is used as a reference to build sync vectors
% raw sync, using theoretical t_btw_trig
trigger_time_begraw = fliplr(Time(idxref2):-t_btw_trig:Time(idxft));
trigger_time_endraw = Time(idxref2):t_btw_trig:Time(idxlt);
% adjusted sync, using actual_t_btw_trig
trigger_time_beg = fliplr(Time(idxref2):-actual_t_btw_trig:Time(idxft));
trigger_time_end = Time(idxref2):actual_t_btw_trig:Time(idxlt);

% to make sure the real first trigger is not missed.
ftmismatch = abs(trigger_time_beg(1) - Time(idxft));
if ftmismatch > 0.5*actual_t_btw_trig 
    trigger_time_beg = [trigger_time_beg(1) - actual_t_btw_trig trigger_time_beg];
    fprintf('First trigger added\n\n');
end

% the sample corresponding to idxref2 is both in "beg" and "end" so we remove it from "beg".
trigger_timeraw = [trigger_time_begraw(1:end-1) trigger_time_endraw];  % raw sync
trigger_time = [trigger_time_beg(1:end-1) trigger_time_end]; % adjusted sync 

check = diff(trigger_timeraw); % check time between triggers (should be constant...)
check2 = diff(trigger_time); % check time between triggers (should be constant...)

% this number should equal the number of acoustic files
Ntriggerraw = length(trigger_timeraw);
Ntrigger = length(trigger_time);

TotalNumberOfFiles = Ntrigger;
fprintf(['The total number of acoustic files should be ' num2str(TotalNumberOfFiles) '.\n\n']);

trigsraw = -0.5*ones(Ntriggerraw,1); % -0.5 to be adjusted depending on the experiment (adjust it to see both signals clearly on figure 2)
trigs = -0.51*ones(Ntrigger,1); 

% plot raw and adjusted sync, as well as the two reference triggers
figure(2)
plot(gca,trigger_timeraw,trigsraw,'*');
plot(gca,trigger_time,trigs,'*');
plot(gca,Time(idxref1),-0.52,'s');
plot(gca,Time(idxref2),-0.52,'s');
set(gca,'FontSize',16);
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
legend('Sync Signal','Raw Sync','Adjusted Sync (using biax recorder as reference)', ...
    'First Reference Trigger','Second Reference Trigger')

% Find peaks in the sync data (will work when recording rate is 1000Hz or higher)
% uncomment to see the experimental peaks...
% [pks,locs] = findpeaks(-Sync(idxft:idxlt),'MinPeakDistance',t_btw_trig*980); 
% % 980 is slightly lower than 1000Hz, the sampling rate used for this run
% locs = locs + idxft - 1;
% plot(Time(locs),Sync(locs),'sm');hold on

figure;
% plot(diff(Time(locs))*1000,'Marker','.');hold on % uncomment to see the
% experimental peaks...

plot(check*1000,'r');hold on
plot(check2*1000,'g');title('Time between consecutive triggers should be constant');
xlabel('Index Number');
ylabel('Time between consecutive triggers (ms)')
legend('Expected','Adjusted (using biax recorder as reference)')
set(gca,'FontSize',14);
fprintf(['Check that the time between consecutive triggers is constant,\n' ...
             'then zoom on the other figure to check that the sync is correct.\n\n'])                          

acN = TotalNumberOfFiles*numWFpfilepCH; % total number of WF per channel

% build acoustic time vector
acTime = (0:acPeriod_adjusted:acPeriod_adjusted*(acN-1))/1e6; % seconds
acTime = acTime + trigger_time(1);                            % seconds

% figure to check that 'acTime' vector is correct (it should overlap 'trigger_time' vector)
% figure; 
% plot(trigger_time,trigger_time,'*');hold on
% plot(acTime(1:numWFpfilepCH:end),acTime(1:numWFpfilepCH:end),'s');
% legend('trigger_time','acTime');
% dcmObj = datacursormode;
% set(dcmObj,'UpdateFcn',@GoodCursor); 

end



