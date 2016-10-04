clear all
close all
clc
% Script to sync acoustic data (verasonics device) with mechanical data

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

acousticrun = 'run2'; % select the acoustic run to sync after adjusting indexes (from figure 1) and the path
switch acousticrun
                
    case 'run1'        
        AcSettingsfile = 'p4640_run1.mat'; % acoustic settings file (located in the current directory)
        WF_path = '/Volumes/DataJRLA/acousticdataPSU/p4640ac/run1/WF_'; % where the WFs are
        % Below, indexes TBD when running the program with step 1        
        idxft = 29714;       % pick up the first trigger of the run 
        idxlt = 2357279;     % pick up the last trigger of the run 
        idxref1 = 1397788;   % pick up a large trigger towards the beginning of the run 
        idxref2 = 2053150;   % pick up a large trigger towards the end of the run                           

    case 'run2'        
        AcSettingsfile = 'p4640_run2.mat'; % acoustic settings file (located in the current directory)
        WF_path = '/Volumes/DataJRLA/acousticdataPSU/p4640ac/run2/WF_'; % where the WFs are
        % Below, indexes TBD when running the program with step 1        
        idxft = 2406190;     % pick up the first trigger of the run 
        idxlt = 3628856;     % pick up the last trigger of the run 
        idxref1 = 2444084;   % pick up a large trigger towards the beginning of the run 
        idxref2 = 3416888;   % pick up a large trigger towards the end of the run           
        
    case 'run3'        
        AcSettingsfile = 'p4640_run3.mat'; % acoustic settings file (located in the current directory)
        WF_path = '/Volumes/DataJRLA/acousticdataPSU/p4640ac/run3/WF_'; % where the WFs are
        % Below, indexes TBD when running the program with step 1        
        idxft = 3659897;       % pick up the first trigger of the run 
        idxlt = 7085194;     % pick up the last trigger of the run 
        idxref1 = 4156544;   % pick up a large trigger towards the beginning of the run 
        idxref2 = 5958790;   % pick up a large trigger towards the end of the run                                 
                
end  

%% plot raw data

% plot Sync signals and other mechanical parameters. 
% From this figure, pick up the first trigger (idxft), the last trigger (idxlt)
% and two reference triggers (idxref1 and idxref2) towards the beginning
% and the end of the run.

FigRaw = figure;
subplot(311);plot(Time,ShearStress);ylabel('Shear Stress (MPa)');
subplot(312);plot(Time,Sync);ylabel('Sync');hold on
subplot(313);plot(Time,NormStress);ylabel('Normal Stress (mm)');
xlabel('Time (s)')
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
fprintf('Pick up indexes from the sync signal (idxft,idxlt,idxref1,idxref2),\n')

% return % uncomment to first pick your indexes

%% sync save data
[acTime,acPeriod,ts,totalnumberoffiles] = SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2);

% save sync data for each run
filenamedata = [runname '_sync_' acousticrun '.mat']; % filename of the mat file
save(filenamedata,'acTime','acPeriod','ts','totalnumberoffiles'); 
              
return
