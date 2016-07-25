clear all
close all
clc
% Main script to analyze acoustic data (verasonics device) and sync it with
% mechanical data.

% Adjust 'Step' to 1, 2, 3, 4, 5 or 6 before running the program
Step = '2';

% Step = '1' to start and choose the run to analyze (e.g. control file X...).
% Step = '2' to check that the sync is correct
% Step = '3' to display WFs at different times during the run (with the option of stacking them)
% Step = '4' to check the reference WF (stack of the first Nstackref waveforms) and
%            pick up the reference time-of-flight and the part of the WF to analyze. 
% Step = '5' to start processing acoustic data and displaying all WFs (slow but
%            useful to check that everything is right).
% Step = '6' to process acoustic data and display only one WFs per acoustic file (faster than 4).

runname = 'p4581';
% read binary file (output of r_file)
[data,outname] = ReadBinBiax(runname);

%% Allocate imported array to column variable names

Time            = data(:,2);
LPDisp          = data(:,3);
ShearStress     = data(:,4);
NormDisp        = data(:,5);
NormStress      = data(:,6);
Sync            = data(:,7);
SampleFreq      = data(:,8);

% ----------------------------------------------------------------------------
% |Column|             Name|             Unit|          Records|
% ----------------------------------------------------------------------------
% |     1|             Time|              sec|          4190772|
% |     2|           LPDisp|              mic|          4190772|
% |     3|      ShearStress|              MPa|          4190772|
% |     4|         NormDisp|           micron|          4190772|
% |     5|       NormStress|              MPa|          4190772|
% |     6|             Sync|              bit|          4190772|
% |     7|        SamplFreq|               Hz|          4190772|
% ----------------------------------------------------------------------------

acousticrun = 'run1'; % select the acoustic run to analyze after adjusting indexes (from figure 1) and the path
switch acousticrun
                
    case 'run1'        
        general_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4581ac';
        run_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4581ac/run1/WF_';
        % Below, indexes TBD when running the program with step 1        
        idxft = 118373;      % pick up the first trigger of the run 
        idxlt = 3460519;     %3461856;   pick up the last trigger of the run 
        idxref1 = 388311;    % pick up a large trigger towards the beginning of the run 
        idxref2 = 3448492;   % pick up a large trigger towards the end of the run 
                        
        % Using step 3, show me WFs at these times
        AtWhichTimes = [2136.85 4737 4741.7]; %3056 3062 vector of times (seconds) at which you would like to see the waveforms
        % AtWhichTimes = [1700:500:3200]; 
        ParamAE.Nsmooth = 400; % moving average Nsmooth samples     
        ParamAE.minPkDist = 0.15e-3; % minimum time interval between 2 consecutive events in seconds
        ParamAE.minPkHeight = 15; % minimum peak height (bit number)      
        ParamAE.threshold = 0.41848e-3; % define as NaN initially. Once setupAEdetection is run, setup
        % this parameter using the threshold output of setupAEdetection.
        % threshold should be in units of seconds (ring-down time). Defining
        % as NaN initially prevents from picking up events everytime you
        % run the program when checking on other parameters.
        
        AnalyzeAllAE = 0; % assign 1 to analyze all AE data, once ParamAE has been setup correctly
        showEvents = 1; % when AnalyzeAllAE is 1, showEvents (1) or not (0)
                
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    No change needed below                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stepoptions = selectstep(Step);

% plot Sync signals and other mechanical parameters. 
% From this figure, pick up the first trigger (idxft), the last trigger (idxlt)
% and two reference triggers (idxref1 and idxref2) towards the beginning
% and the end of the run.

FigRaw = figure;
subplot(311);plot(Time,ShearStress);ylabel('Shear Stress (MPa)');
subplot(312);plot(Time,Sync);ylabel('Sync');hold on
subplot(313);plot(Time,NormStress);ylabel('Normal Stress (MPa)');
xlabel('Time (s)')
dcmObj = datacursormode;set(dcmObj,'UpdateFcn',@GoodCursor);
if stepoptions.CHOOSE_RUN
    fprintf(['Pick up indexes from the sync signal (idxft,idxlt,idxref1,idxref2),\n then run the program again with step 2.\n\n'])
    break
end

filenamedata = [acousticrun 'ac.mat']; % filename of the mat file
AcSettingsfile = [runname '_AE.mat'];

% sync data
[acTime,acRate_adjusted,ts_adjusted,totalnumberoffiles] = SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2);
if stepoptions.CHECK_SYNC
    break
end

if AnalyzeAllAE == 1
    % analyze all data
    [EvAmp,EvTime] = AE(AcSettingsfile,run_ac_path,ts_adjusted,acTime,totalnumberoffiles,ParamAE,showEvents);
else
    % prepare AE parameters
    RD = setupAEdetection(AcSettingsfile,run_ac_path,ts_adjusted,acTime,AtWhichTimes,ParamAE);
    return
end






%% save data
save(filenamedata,'acTime','EvAmp','EvTime');                                                           
              
