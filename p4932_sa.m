clear all
close all
clc
% Script to sync acoustic data (verasonics device) with mechanical data

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


acousticrun = 'run1'; % select the acoustic run to sync after adjusting indexes (from figure 1) and the path
switch acousticrun                                

    case 'run1'        
        AcSettingsfile = 'p4932_run1.mat'; % acoustic settings file (located in the current directory)
        WF_path = '/Volumes/JAQUEZ/p4932ac/run1/WF_'; % where the WFs are
%         totalnumberoffiles = 2302;
        % Below, indexes TBD when running the program with step 1        
        idxft = 219826;     % pick up the first trigger of the run 
        idxlt = 551312;     %61331 pick up the last trigger of the run 
        idxref1 = 305361;   % pick up a large trigger towards the beginning of the run 
        idxref2 = 464769;   % pick up a large trigger towards the end of the run   
    
    case 'run2'        
        AcSettingsfile = 'p4932_run2.mat'; % acoustic settings file (located in the current directory)
        WF_path = '/Volumes/JAQUEZ/p4932ac/run2/WF_'; % where the WFs are
        totalnumberoffiles = 1008;                                                  
end  

%% plot raw data

% plot Sync signals and other mechanical parameters. 
% From this figure, pick up the first trigger (idxft), the last trigger (idxlt)
% and two reference triggers (idxref1 and idxref2) towards the beginning
% and the end of the run.

FigRaw = figure;
ax1 = subplot(311);plot(Time,ShearStress);ylabel('Shear Stress (MPa)');
ax2 = subplot(312);plot(Time,Sync);ylabel('Sync');hold on % Sync
ax3 = subplot(313);plot(Time,NormStress);ylabel('Normal Stress (MPa)');
xlabel('Time (s)')
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
fprintf('Pick up indexes from the sync signal (idxft,idxlt,idxref1,idxref2),\n')

linkaxes([ax1,ax2,ax3],'x');

% return
% uncomment to first pick your indexes

%% sync save data
[acTime,acPeriod,ts,totalnumberoffiles] = SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2);

% save sync data for each run
filenamedata = [runname '_sync_' acousticrun '.mat']; % filename of the mat file
save(filenamedata,'acTime','acPeriod','ts','totalnumberoffiles'); 
              
return