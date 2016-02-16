clear all
close all
clc
% Main script to analyze acoustic data (verasonics device) and sync it with
% mechanical data.

% Adjust 'Step' to 1, 2, 3, 4, 5 or 6 before running the program
Step = '3';

% Step = '1' to start and choose the run to analyze (e.g. control file X...).
% Step = '2' to check that the sync is correct
% Step = '3' to display WFs at different times during the run (with the option of stacking them)
% Step = '4' to check the reference WF (stack of the first 50 waveforms) and
%            pick up the reference time-of-flight and the part of the WF to analyze. 
% Step = '5' to start processing acoustic data and displaying all WFs (slow but
%            useful to check that everything is right).
% Step = '6' to process acoustic data and display only one WFs per acoustic file (faster than 4).

% experiment #
runname = 'p4530';
% read binary file (output of r_file)
[data,outname] = ReadBinBiax(runname);

% adjust name and number depending on the experiment
SampleNum    = data(:,1);
LPDisp       = data(:,2);
ShearStress  = data(:,3);
NormDisp     = data(:,4);
NormStress   = data(:,5);
Time         = data(:,6);
ec_dispn     = data(:,7);
layer_thick  = data(:,8);
ec_dispv     = data(:,9);
SampleFreq   = data(:,10);
mu           = data(:,11);
Shear_Strain = data(:,12);
Sync         = data(:,13);

% ----------------------------------------------------------------------------
% |Column|             Name|             Unit|          Records|
% ----------------------------------------------------------------------------
% |     1|          LP_Disp|              mic|           271431|
% |     2|       Shr_stress|              MPa|           271431|
% |     3|         nor_disp|           micron|           271431|
% |     4|       Nor_stress|              MPa|           271431|
% |     5|             Time|              sec|           271431|
% |     6|         ec_dispn|           micron|           271431|
% |     7|      layer_thick|           micron|           271431|
% |     8|         ec_dispv|              mic|           271431|
% |     9|        Samp_Freq|               Hz|           271431|
% |    10|               mu|                .|           271431|
% |    11|     Shear_Strain|                .|           271431|
% |    12|             Sync|              bit|           271431|
% ----------------------------------------------------------------------------

acousticrun = 'run1_P'; % select the acoustic run to analyze after adjusting indexes (from figure 1) and the path
switch acousticrun
                
    case 'run1_P'
        general_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4530';
        run_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4530/run1/WF_';
        % Below, indexes TBD when running the program with step 1
        idxft = 31967; % pick up the first trigger of the run 
        idxlt = 259798; % pick up the last trigger of the run 
        idxref1 = 35007; % pick up a large trigger towards the beginning of the run 
        idxref2 = 249568; % pick up a large trigger towards the end of the run 
        
        % Below, indexes TBD when running the program with step 4
        idxBeg = 442;       % choose the beginning of the WF used for analysis
        idxEnd = 707;       % choose the end of the WF used for analysis
        idx_TOF_0 = 442;    % idx corresponding to the arrival of the WFref 
        
        % number of waveforms to stack (either to display or when analyzing)
        NtoStack = 500;
        
    case 'run1_S'
        general_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4530';
        run_ac_path = '/Volumes/DataJRLA/acousticdataPSU/p4530/run1/WF_';
        % Below, indexes TBD when running the program with step 1
        idxft = 31967; % pick up the first trigger of the run 
        idxlt = 259798; % pick up the last trigger of the run 
        idxref1 = 35007; % pick up a large trigger towards the beginning of the run 
        idxref2 = 249568; % pick up a large trigger towards the end of the run                 
        
        % Below, indexes TBD when running the program with step 4
        idxBeg = 764;      % choose the beginning of the WF used for analysis
        idxEnd = 892;      % choose the end of the WF used for analysis
        idx_TOF_0 = 764;   % idx corresponding to the arrival of the WFref   

        % number of waveforms to stack (either to display or when analyzing)
        NtoStack = 20; 

end  

% Using step 3, show me WFs at these times
% AtWhichTimes = [1700 2000 2500 3000]; % vector of times (seconds) at which you would like to see the waveforms
AtWhichTimes = [1700:500:3200];
% choose 'absoluteref' to cross-correlate all WFs with the NtoStackref first ones. 
% or choose 'relativeref' to cross-correlate each WF with the previous one.
ref = 'absoluteref'; %'absoluteref';
% when using 'absoluteref', NtoStackref is the number of WFs used to build a
% reference WF. When using 'relativeref', this parameter is set equal to 'NtoStack'
NtoStackref = 100; 

% threshold is used to ignore noisy waveforms when using 'relativeref'.
% TimeShift values are kept only when MaxInter is above this threshold.
% When 'absoluteref' is used, threshold is ignored.
threshold = 0.94;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    No change needed below                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stepoptions = selectstep(Step);

% plot Sync signals and other mechanical parameters. 
% From this figure, pick up the first trigger (idxft), the last trigger (idxlt)
% and two reference triggers (idxref1 and idxref2) towards the beginning
% and the end of the run.

FigRaw = figure;
subplot(411);semilogy(Time,SampleFreq);ylabel('Sampl. Freq.');
subplot(412);plot(Time,Sync);ylabel('Sync');hold on
subplot(413);plot(Time,LPDisp);ylabel('LP disp (microns)');
subplot(414);plot(Time,ShearStress);ylabel('Shear Stress (MPa)');
xlabel('Time (s)')
dcmObj = datacursormode;
set(dcmObj,'UpdateFcn',@GoodCursor);
fprintf(['Pick up indexes from the sync signal (idxft,idxlt,idxref1,idxref2),\n then run the program again with step 2.\n'])
if stepoptions.CHOOSE_RUN
    break
end

filenamedata = [acousticrun 'ac.mat']; % filename of the mat file
AcSettingsfile = [general_ac_path '/' runname '.mat'];

% sync data
[acTime,acRate_adjusted,ts_adjusted,totalnumberoffiles] = SyncAcData(AcSettingsfile,Time,Sync,idxft,idxlt,idxref1,idxref2);
if stepoptions.CHECK_SYNC
    break
end

% show WFs at different times
ShowMeWFs(AcSettingsfile,run_ac_path,ts_adjusted,acTime,AtWhichTimes,NtoStack);
if stepoptions.SHOW_WFbefAnalysis
    break
end

% process acoustic data (Time Shift, RmsAmp and Max Intercorrelation)
[MaxInter,TimeShift,RmsAmp,Amp,TOF_0,RmsAmpRef,AmpRef,fullWFref] = ...
    ProcessAc_newtest(AcSettingsfile,run_ac_path,ts_adjusted,totalnumberoffiles,...
    idxBeg,idxEnd,idx_TOF_0,stepoptions,ref,NtoStackref,NtoStack,threshold); %_newtest
if stepoptions.SHOW_WFref
    break
end

if strcmp(ref,'relativeref')
    idxNaN = isnan(TimeShift); % find NaNs
    TimeShift(idxNaN) = 0; % 0 instead of NaN to be able to use cumsum below
    
    % sum up all the delays and put back NaN into the final vector
    TimeShift = cumsum(TimeShift);
    TimeShift(idxNaN) = NaN;
end

figure;
subplot(311);plot(TimeShift);
subplot(312);plot(MaxInter);
subplot(313);plot(RmsAmp);

%% save data
save(filenamedata,...
    'acTime','RmsAmp','Amp','AmpRef','MaxInter','TimeShift',...        
    'TOF_0','RmsAmpRef',...                             
    'fullWFref','idxBeg','idxEnd','NtoStackref','NtoStack','ref','threshold');                                                           
              
return
