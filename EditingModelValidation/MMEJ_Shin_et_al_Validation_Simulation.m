clear all
%% Villiger et al validation simulation

% Load simulation
sbioloadproject('SimBiologyProjectFiles\FAH_Editing_Project_Shin_2018.sbproj')
csObj = getconfigset(m1);
DoseObj = getdose(m1);

%% Edit the simbiology project file to match the experimental circumstances
% NTBC withdrawn on day 7
m1.Events(1).Trigger = 'time>7';
% csObj(1).StopTime = 37;
m1.Species(8).Value = m1.Species(8).Value / 2; % Reduce DNA dose to match the 30 ug dosing  
%% Get the list of editing rate constants from the DNA editing models
% Load Validation Data structure and calculate the mouse RNP dose
load('ValidationdEdtData/ValidationStructureFile_dEdt_Norm.mat')
k_list = [ValidationModel(1).BootStrappedk ValidationModel(2).BootStrappedk];
k_total_edit = [min(k_list) mean(k_list) max(k_list)];

csObj = getconfigset(m1);
DoseObj = getdose(m1);
csObj(1).StopTime = 37;

for i=1:3
    [Time, Cells(i).Cells, RNP(i).RNP, dE_dt_Norm(i).dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(i));
end

for i=1:3
    PreciseEditRange(i) = Cells(i).Cells.Precise(end) / sum([Cells(i).Cells.Precise(end) Cells(i).Cells.Imprecise(end) Cells(i).Cells.Unedited(end)]);
end

csObj(1).StopTime = 35;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
set(gca,'YLim', [0 .20]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Shin_Validation_TimeCourse.fig');
figure;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, dEdt_Norm);
set(gca,'YLim', [0 .120]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Villiger_Shin_dEdt_Norm.fig');

%% Extracted Villiger Data for plotting


PublishedData_Mean = 5.18;
PublishedData_Top = 5.18 + 1.92;
PublishedData_Bottom = 5.18 - 1.92;

