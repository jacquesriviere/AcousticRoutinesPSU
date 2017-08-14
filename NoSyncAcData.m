function [acTime,acPeriod,ts] = NoSyncAcData(Ac_path,TotalNumberOfFiles)
% NoSyncAcData uses acoustic settings (mat file from Verasonics) and the
% number of files provided by user to create a time vector starting at t =
% 0 and having one datapoint per waveform. This function is used when
% acoustic data is not synchronized with any other data (e.g. biax data).
% It replaces SyncAcData.

% The main output is a time vector for acoustic data. It also outputs an
% adjusted acoustic sampling rate and an adjusted acoustic pulsing rate.

% INPUTS
% Ac_path is the path where pXXXX.mat (containing acoustic settings) can be
% loaded
% TotalNumberOfFiles is the total number of acoustic files to be analyzed.

% OUTPUTS
% acTime is the time vector for acoustic data obtained after synchronization
% acRate is the rate at which acoustic pulses are sent 
% ts is the sampling period in microsec (i.e. 1/fs)

% acoustic parameters
acSettings = load(Ac_path);                     % load acoustic settings
acPeriod = acSettings.SeqControl(1).argument;   % time btw pulses in microsec (SeqControl(1).argument = timeBetweenpulses)
numSFpfile = acSettings.numFrames/2;            % number of superframes per file
numWFpSFpCH = acSettings.numAcqs;               % number of WF per superframe and per channel
numWFpfilepCH = numSFpfile*numWFpSFpCH;         % number of WF per file and per channel
fs = acSettings.samplingFreq;                   % acoustic sampling rate in MHz
clear acSettings

ts = 1/fs; % sampling time in microsec                        
         
acN = TotalNumberOfFiles*numWFpfilepCH; % total number of WF per channel

% build acoustic time vector
acTime = (0:acPeriod:acPeriod*(acN-1))/1e6; % seconds

fprintf(['The acoustic time vector goes from 0 to ' num2str(acTime(end)) ' s.\n']);

end
