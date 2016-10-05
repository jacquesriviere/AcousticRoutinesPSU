clear all
close all
clc

%% load acoustic data (output from pXXXX_a.m)

% abs
run = 'Results_p4640_run2_10468.0174s-11034.9828s_Stack10WFs_absref';
load([run '.mat']);

TimeShiftAbs = TimeShift;
RmsAmpAbs = RmsAmp;
MaxInterAbs = MaxInter;

% rel
run = 'Results_p4640_run2_10468.0174s-11034.9828s_Stack10WFs_relref_Th0.95';
load([run '.mat']);

TimeShiftRel = TimeShift;
RmsAmpRel = RmsAmp;
MaxInterRel = MaxInter;

% mix
run = 'Results_p4640_run2_10468.0174s-11034.9828s_Stack10WFs_mixref_Th0.95';
load([run '.mat']);


%% load mechanical data (only in the range where the acoustic run was made)
runname = 'p4640';
% read binary file (output of r_file)
[data,outname] = ReadBinBiax(runname);

Time = data(:,2); % load the entire time vector to figure out which part we need to keep

% idx numbers corresponding to the beginning and end of the acoustic run
FirstIdxAc = find(Time > LocalAcTime(1,1),1,'first'); 
LastIdxAc = find(Time < LocalAcTime(end,end),1,'last');
idxAc = FirstIdxAc:LastIdxAc;

Time            = Time(idxAc); % keep only the useful part of the vectors
LPDisp          = data(idxAc,3);
ShearStress     = data(idxAc,4);
NormDisp        = data(idxAc,5);
NormStress      = data(idxAc,6);
Pc_disp         = data(idxAc,7);
Pc              = data(idxAc,8);
Ppa_disp        = data(idxAc,9);
Ppa             = data(idxAc,10);
Ppb_disp        = data(idxAc,11);
Ppb             = data(idxAc,12);
IntDisp         = data(idxAc,13);
Sync            = data(idxAc,14);
SamplFreq       = data(idxAc,15);
effNorStress    = data(idxAc,16);
Qa              = data(idxAc,17);
Qb              = data(idxAc,18);
Qdiff           = data(idxAc,19);
Qavg            = data(idxAc,20);
PpDiff          = data(idxAc,21);
Perm            = data(idxAc,22);
clear data

%% Display

R = 2; % receiver 2
T = 3; % transmitter 3

% offset mechanical and acoustic time vectors to start the plot at 0.
MechTime = Time - LocalAcTime(1,T);
AcTime = LocalAcTime(:,T) - LocalAcTime(1,T);

% empty plot with 2 y-axes
empty = NaN(10,1);
figure; 
[ax, p1, p2]=plotyy(empty,empty,empty,empty); 
delete(p1);delete(p2);

hold(ax(1),'on')
scatter(ax(1),AcTime,squeeze(TimeShiftAbs(:,R,T)),[],squeeze(MaxInterAbs(:,R,T)),'.');
set(gca,'CLim',[0.5 1]);
colorbar('northoutside');
scatter(ax(1),AcTime,squeeze(TimeShiftRel(:,R,T)),[],squeeze(MaxInterRel(:,R,T)),'.');
scatter(ax(1),AcTime,squeeze(TimeShift(:,R,T)),[],squeeze(MaxInter(:,R,T)),'.');
set(ax(1),'XLimMode','auto','YLimMode','auto');
set(ax(1),'XTickMode','auto','YTickMode','auto');
set(ax(1),'FontSize',16);
xlabel('Time (s)','Interpreter','Latex');
ylabel('Time Shift ($\mu$s)','Interpreter','Latex');

hold(ax(2),'on')
plot(ax(2),MechTime,ShearStress,'LineWidth',2) 
set(ax(2),'XLim',get(ax(1),'XLim'),'YLimMode','auto'); 
set(ax(2),'XTick',get(ax(1),'XTick'),'YTickMode','auto');
set(ax(2),'FontSize',16);
ylabel(ax(2),'Shear Stress (MPa)','Interpreter','Latex');
