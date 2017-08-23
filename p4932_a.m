clear all
close all
clc
% Script to anaylize acoustic data (verasonics device), once the sync has
% been performed.

% Load mechanical data (filename looks like 'pXXXX_data_bin')
runname = 'p4932';

% read binary file (output of r_file)
[data,outname] = ReadBinBiax(runname);

LPDisp          = data(:,2);
ShearStress     = data(:,3);
NormDisp        = data(:,4);
NormStress      = data(:,5);
Time            = data(:,6);
Sync            = data(:,7);
Samp_Freq       = data(:,8);
mu              = data(:,9);

% ----------------------------------------------------------------------------
% |Column|             Name|             Unit|          Records|
% ----------------------------------------------------------------------------
% |     1|          LP_Disp|              mic|           733509|
% |     2|       Shr_stress|              MPa|           733509|
% |     3|         nor_disp|           micron|           733509|
% |     4|       Nor_stress|              MPa|           733509|
% |     5|             Time|              sec|           733509|
% |     6|             sync|              bit|           733509|
% |     7|        Samp_Freq|               Hz|           733509|
% |     8|               mu|                .|           733509|
% ----------------------------------------------------------------------------

acousticrun = 'run1'; % select the acoustic run to analyze after adjusting indexes (from figure 1) and the path
switch acousticrun  

    case 'run1'        
        AcSettingsfile = 'p4932_run1.mat'; % acoustic settings file (located in the current directory)
        AcSyncFile = 'p4932_sync_run1.mat'; % sync file                                 
        WF_path = '/Volumes/JAQUEZ/p4932ac/run1/WF_'; % where the WFs are     
        % Portion of the WF to analyze
        idxBeg = 376;       % choose the beginning of the WF used for analysis
        idxEnd = 587;       % choose the end of the WF used for analysis                
        % Display waveforms sent by transmitter WhichTrans
        WhichTrans = 1;           
        % Time Range
        TimeRange = [4500 6000]; % in seconds
        AtWhichTimes = linspace(4500,6000,6); % first set of oscillations (post fracture)    
    
    case 'run2'        
        AcSettingsfile = 'p4932_run2.mat'; % acoustic settings file (located in the current directory)
        AcSyncFile = 'p4932_nosync_run2.mat'; % sync file                                 
        WF_path = '/Volumes/JAQUEZ/p4932ac/run2/WF_'; % where the WFs are     
        % Portion of the WF to analyze
        idxBeg = 462;       % choose the beginning of the WF used for analysis
        idxEnd = 608;       % choose the end of the WF used for analysis                
        % Display waveforms sent by transmitter WhichTrans
        WhichTrans = 3;           
        % Time Range
%         TimeRange = [500 700];
        AtWhichTimes = linspace(100,1000,6); % first set of oscillations (post fracture)
%         AtWhichTimes = 7207.2; % first set of oscillations (post fracture)
        
end  

% number of waveforms to stack (either to display or when analyzing)
NtoStack = 10;

% analyze acoustic data over that time range only. If not defined,
% the whole run is analyzed

%         TimeRange = [10468 11035]; % time interval in seconds

displayoptions = 1; % choose 0 to display all waveforms or 1 to display one set of waveforms over 100

% Show WFs at these times
% vector of times (seconds) at which you would like to see the waveforms

% AtWhichTimes = linspace(13500,15000,10); % around first slip

ref = 'absref'; %'absref', 'relref' or 'mixref';

% offset waveforms by Offset when multiple channels are used
Offset1 = 10000;
Offset2 = 10000;

% used for 'relref' or 'mixref'
threshold = 0.95;

%% plot mechanical data of interest within the considered run

load(AcSyncFile);
% Find sample number corresponding to the beginning and end of the acoustic run
FirstIdxAc = find(Time > acTime(1),1,'first'); 
LastIdxAc = find(Time > acTime(end),1,'first');
idxAc = FirstIdxAc:LastIdxAc;

FigRaw = figure;
ax1 = subplot(311);plot(Time(idxAc),ShearStress(idxAc));ylabel('Shear Stress (MPa)');
ax2 = subplot(312);plot(Time(idxAc),LPDisp(idxAc)/1000);ylabel('LP Disp (mm)');hold on
ax3 = subplot(313);plot(Time(idxAc),NormStress(idxAc));ylabel('Normal Stress (MPa)');
xlabel('Time (s)')
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);

linkaxes([ax1,ax2,ax3],'x');


%% show WFs at different times
ShowMeWFs(WF_path,AcSettingsfile,AcSyncFile,AtWhichTimes,NtoStack,Offset1,WhichTrans);

%% process acoustic data (Time Shift, RmsAmp and Max Intercorrelation)
[MaxInter,TimeShift,RmsAmp,Amp,RmsAmpRef,AmpRef,fullWFref,LocalAcTime] = ...
    ProcessAc_Tomo(WF_path,AcSettingsfile,AcSyncFile,...
    idxBeg,idxEnd,ref,NtoStack,threshold,Offset2,displayoptions,TimeRange); %,TimeRange
%% save data

% filename of the resulting mat file with explicit name based on chosen paramters
if strcmp(ref,'relref') || strcmp(ref,'mixref')
   ref = [ref '_Th' num2str(threshold)]; % add threshold value to name when using relative reference
end
try % if TimeRange is defined above
    filenamedata = ['Results_' runname '_' acousticrun '_' num2str(LocalAcTime(1,1)) 's-' num2str(LocalAcTime(end,end)) 's_Stack' num2str(NtoStack) 'WFs_' ref '.mat'];
    save(filenamedata,...
    'LocalAcTime','RmsAmp','Amp','AmpRef','MaxInter','TimeShift',...        
    'RmsAmpRef','fullWFref','idxBeg','idxEnd','NtoStack','ref','threshold','TimeRange');

catch % if TimeRange is not defined (full run analyzed)
    filenamedata = ['Results_' runname '_' acousticrun '_fullrun_Stack' num2str(NtoStack) 'WFs_' ref '.mat'];
    TimeRange = [acTime(1) acTime(end)];
    save(filenamedata,...
    'LocalAcTime','RmsAmp','Amp','AmpRef','MaxInter','TimeShift',...        
    'RmsAmpRef','fullWFref','idxBeg','idxEnd','NtoStack','ref','threshold','TimeRange');
end
    
                                                           
              
return
