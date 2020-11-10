clear all
%% Villiger et al validation simulation

% Load simulation
sbioloadproject('SimBiologyProjectFiles\PAH_Editing_Villiger_et_al_2018.sbproj')

%% Editing the AAV simulation scheme to suit Villiger
% % Change the DNA dose
% m1.Species(8).InitialAmount = 2.843e9 / 1e12 * 5e11; %Number is based on 1e11 vector genomes yielding 2.843e9 genome copies of editor DNA
% 
% % Change the NTBC toggle to zero, essentially avoiding any selection
% m1.Events(1).EventFcns = 'NTBC_Removal_Toggle = 0';
% 
% % Code in the specificity by altering the precise imprecise ratio
% FractionPrecise = 40/(22+40); % Precision based on data from Figure 1d
% k_Total_Editing = 1.6857e-13; % high 1.6857e-13 low 2.5498e-14
% m1.Parameters(8).Value = FractionPrecise * k_Total_Editing;
% m1.Parameters(9).Value = (1 - FractionPrecise) * k_Total_Editing;
% 
% % Change simulation time to 4 weeks
% csObj = getconfigset(m1);
% csObj(1).StopTime = 3 * 7;
% 
% sbiosaveproject('SimBiologyProjectFiles\PAH_Editing_Villiger_2018.sbproj', 'm1');
%% Get the list of editing rate constants from the DNA editing models
% Load Validation Data structure and calculate the mouse RNP dose
load('ValidationdEdtData/ValidationStructureFile_dEdt_Norm.mat')

% Calculate average and standard deviation of dEdt_Norm
for i=1:4
    dEdt_Data(i).Mean = mean(ValidationModel(i).NormalizeddEdt);
    dEdt_Data(i).Std  = std(ValidationModel(i).NormalizeddEdt);
end

Aggregate_dEdt_Mean = mean([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
Aggregate_dEdt_Std = std([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
Plasmid_dEdt_Mean = mean([dEdt_Data(1).Mean, dEdt_Data(2).Mean]);
Plasmid_dEdt_Std = std([dEdt_Data(1).Mean, dEdt_Data(2).Mean]);

dEdt_Range = [(Aggregate_dEdt_Mean-Aggregate_dEdt_Std) Aggregate_dEdt_Mean (Aggregate_dEdt_Mean+Aggregate_dEdt_Std)];

% Iterate to find k
for i=1:3
    [k_total_edit(i) k_Modifier k_History] = Generic_k_Finder_AvgdEdt(m1, dEdt_Range(i));
end

csObj = getconfigset(m1);
DoseObj = getdose(m1);
csObj(1).StopTime = 21;

for i=1:3
    [Time, Cells(i).Cells, RNP(i).RNP, dE_dt_Norm(i).dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(i));
end

for i=1:3
    PreciseEditRange(i) = Cells(i).Cells.Precise(end) / sum([Cells(i).Cells.Precise(end) Cells(i).Cells.Imprecise(end) Cells(i).Cells.Unedited(end)]);
end

csObj(1).StopTime = 35;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
set(gca,'YLim', [0 .20]);
savefig('FigurePlots\Villiger_Validation_TimeCourse.fig');
figure;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, dEdt_Norm);
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, dEdt_Norm);
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, dEdt_Norm);
set(gca,'YLim', [0 .120]);
savefig('FigurePlots\Villiger_Validation_dEdt_Norm.fig');

%% Extracted Villiger Data for plotting
DigitizedData(1) = 8.035714285714292;
DigitizedData(2) = 14.28571428571429;
DigitizedData(3) = 17.410714285714292;
DigitizedData(4) = 20.535714285714292;

PublishedData_Mean = DigitizedData(3) - DigitizedData(1);
PublishedData_Top = DigitizedData(4) - DigitizedData(1);
PublishedData_Bottom = DigitizedData(2) - DigitizedData(1);

