clear all
close all
clc
% Script to anaylize acoustic data (verasonics device), once the sync has been performed

% Load mechanical data (filename looks like 'pXXXX_data_bin')
runname = 'p4640';
% read binary file (output of r_file)
[data,outname] = ReadBinBiax(runname);

%% Allocate imported array to column variable names

Time            = data(:,2);
LPDisp          = data(:,3);
ShearStress     = data(:,4);
NormDisp        = data(:,5);
NormStress      = data(:,6);
Pc_disp         = data(:,7);
Pc              = data(:,8);
Ppa_disp        = data(:,9);
Ppa             = data(:,10);
Ppb_disp        = data(:,11);
Ppb             = data(:,12);
IntDisp         = data(:,13);
Sync            = data(:,14);
SamplFreq       = data(:,15);
effNorStress    = data(:,16);
Qa              = data(:,17);
Qb              = data(:,18);
Qdiff           = data(:,19);
Qavg            = data(:,20);
PpDiff          = data(:,21);
perm            = data(:,22);
clear data
% return
% ----------------------------------------------------------------------------
% |Column|             Name|             Unit|          Records|
% ----------------------------------------------------------------------------
% |     1|             Time|              sec|          7125873|
% |     2|           LPDisp|              mic|          7125873|
% |     3|      ShearStress|              MPa|          7125873|
% |     4|         NormDisp|           micron|          7125873|
% |     5|       NormStress|              MPa|          7125873|
% |     6|          Pc_disp|           micron|          7125873|
% |     7|               Pc|              MPa|          7125873|
% |     8|         Ppa_disp|           micron|          7125873|
% |     9|              Ppa|              MPa|          7125873|
% |    10|         Ppb_disp|           micron|          7125873|
% |    11|              Ppb|              MPa|          7125873|
% |    12|          IntDisp|           micron|          7125873|
% |    13|             Sync|              bit|          7125873|
% |    14|        SamplFreq|               Hz|          7125873|
% |    15|     effNorStress|              MPa|          7125873|
% |    16|               Qa|            m^3/s|          7125873|
% |    17|               Qb|            m^3/s|          7125873|
% |    18|            Qdiff|            m^3/s|          7125873|
% |    19|             Qavg|            m^3/s|          7125873|
% |    20|           PpDiff|               Pa|          7125873|
% |    21|             perm|              m^2|          7125873|
% ----------------------------------------------------------------------------

acousticrun = 'run2'; % select the acoustic run to analyze after adjusting indexes (from figure 1) and the path
switch acousticrun                        
                
    case 'run2'        
        AcSettingsfile = 'p4640_run2.mat'; % acoustic settings file 
        AcSyncFile = 'p4640_sync_run2.mat'; % sync file 
        WF_path = 'p4640ac/run2/WF_'; % where the WFs are
        
        % Portion of the WF to analyze
        idxBeg = 87;       % choose the beginning of the WF used for analysis
        idxEnd = 1024;       % choose the end of the WF used for analysis        
        
        % number of waveforms to stack (either to display or when analyzing)
        NtoStack = 10;      
        
        % analyze acoustic data over that time range only. If not defined,
        % the whole run is analyzed
        TimeRange = [10750 10800];
%         TimeRange = [10468 11035]; % time interval in seconds
        
        displayoptions = 1; % choose 0 to display all waveforms or 1 to display one set of waveforms over 100
        
        % Show WFs at these times
        % vector of times (seconds) at which you would like to see the waveforms
        AtWhichTimes = [9968 10487 10653 10725 10798]; 

        % Display waveforms sent by transmitter WhichTrans
        WhichTrans = 3;
        
        ref = 'mixref'; %'absref', 'relref' or 'mixref';                             
        
end  

% offset waveforms by Offset when multiple channels are used
Offset1 = 5000;
Offset2 = 5000;

% used for 'relref' or 'mixref'
threshold = 0.95;

%% plot mechanical data of interest within the considered run

load(AcSyncFile);
% Find sample number corresponding to the beginning and end of the acoustic run
FirstIdxAc = find(Time > acTime(1),1,'first'); 
LastIdxAc = find(Time > acTime(end),1,'first');
idxAc = FirstIdxAc:LastIdxAc;

FigRaw = figure;
subplot(311);plot(Time(idxAc),ShearStress(idxAc));ylabel('Shear Stress (MPa)');
subplot(312);plot(Time(idxAc),LPDisp(idxAc)/1000);ylabel('LP Disp (mm)');hold on
subplot(313);plot(Time(idxAc),NormStress(idxAc));ylabel('Normal Stress (MPa)');
xlabel('Time (s)')
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);

% return

%% show WFs at different times
ShowMeWFs(WF_path,AcSettingsfile,AcSyncFile,AtWhichTimes,NtoStack,Offset1,WhichTrans);

% return

%% process acoustic data (Time Shift, RmsAmp and Max Intercorrelation)
[MaxInter,TimeShift,RmsAmp,Amp,RmsAmpRef,AmpRef,fullWFref,LocalAcTime] = ...
    ProcessAc_Tomo(WF_path,AcSettingsfile,AcSyncFile,...
    idxBeg,idxEnd,ref,NtoStack,threshold,Offset2,displayoptions,TimeRange); %
%% save data

% filename of the resulting mat file with explicit name based on chosen paramters
if strcmp(ref,'relref') || strcmp(ref,'mixref')
   ref = [ref '_Th' num2str(threshold)]; % add threshold value to name when using relative reference
end
try % if TimeRange is defined above
    filenamedata = ['Results_' runname '_' acousticrun '_' num2str(LocalAcTime(1,1)) 's-' num2str(LocalAcTime(end,end)) 's_Stack' num2str(NtoStack) 'WFs_' ref '.mat'];
catch % if TimeRange is not defined (full run analyzed)
    filenamedata = ['Results_' runname '_' acousticrun '_fullrun_Stack' num2str(NtoStack) 'WFs_' ref '.mat'];
    TimeRange = [acTime(1) acTime(end)];
end
    
save(filenamedata,...
    'LocalAcTime','RmsAmp','Amp','AmpRef','MaxInter','TimeShift',...        
    'RmsAmpRef','fullWFref','idxBeg','idxEnd','NtoStack','ref','threshold','TimeRange');                                                           
              
return
