clear all
%% Validation Script for ICAM-2 editing in Sago et al, Dahlman lab
% Load mRNA editing file
sbioloadproject('SimBiologyProjectFiles\ICAM_2_Editing_Sago_2018.sbproj')

% Assuming mice having weight that is the average of 5 and 12 week
ms_wt = (17.8 + 28.9)/2; % g
Cas9_mRNA_Dose = 1/2 * ms_wt/1000 * 1e-3; %g
Cas9_mRNA_MW = (1368 * 320.5) + 159; % g/mol
Cas9_mRNA_mol = Cas9_mRNA_Dose / Cas9_mRNA_MW;
Cas9_mRNA_molecule = Cas9_mRNA_mol * 6.023e23;

% Provide 2 doses of mRNA using dose rather than initial amount
m1.Species(9).InitialAmount = 0;

DoseObj = getdose(m1,'repeat');
DoseObj(1).Amount = Cas9_mRNA_molecule/100;

% Set stop time to 7 days
csObj = getconfigset(m1);
csObj(1).StopTime = 7;

% Set precise editing rate constant to zero
m1.Parameters(8).Value = 0;

% Import k values from the mRNA simulation
load('ValidationdEdtData\ValidationStructureFile_dEdt_Norm.mat')
% Use the dEdt data

for i=1:4
    dEdt_Data(i).Mean = mean(ValidationModel(i).NormalizeddEdt);
    dEdt_Data(i).Std  = std(ValidationModel(i).NormalizeddEdt);
end

Aggregate_dEdt_Mean = mean([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
Aggregate_dEdt_Std = std([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);

dEdt_Range = [Aggregate_dEdt_Mean - Aggregate_dEdt_Std Aggregate_dEdt_Mean Aggregate_dEdt_Mean + Aggregate_dEdt_Std];

% Cannot use generic functions due to dose object being used, using custom
% functions defined at the end of the file
k_imprecise_edit = m1.Parameters(9).Value;

[Time SimData Names] = sbiosimulate(m1, csObj(1), DoseObj(1));

for i=1:3
    Desired_dEdt = dEdt_Range(i);
    for j=1:100
        m1.Parameters(9).Value = k_imprecise_edit;
        [Time SimData Names] = sbiosimulate(m1, csObj(1), DoseObj(1));
        dEdt_Avg = dEdtCalculator(k_imprecise_edit, Time, SimData, Names);
        k_modifier = Desired_dEdt/dEdt_Avg;
        k_imprecise_edit = k_imprecise_edit * k_modifier;
        [i, j];
    end
    k_total_edit(i) = k_imprecise_edit;
end

%k_total_edit = [2.92053469566353e-15,1.25524009335604e-14,2.34626882316588e-14];

for i=1:3
    m1.Parameters(9).Value = k_total_edit(i);
    [Time, SimData, Names] = sbiosimulate(m1,csObj(1),DoseObj(1));
    N_Datapoints = length(Time);
    TotalCells(i) = SimData(N_Datapoints,1) + SimData(N_Datapoints,2) + SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditedCells(i) = SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditFraction(i) = EditedCells(i)/TotalCells(i);
end
EditFractionForPlot = EditFraction(1) + rand(10,1)*(EditFraction(3) - EditFraction(1));
save('ValidationdEdtData\ICAM_2_Editing_Predictions.mat','EditFraction');

%% Graphing this 
csObj = getconfigset(m1);

csObj(1).StopTime = 35;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, Cells.Imprecise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, Cells.Imprecise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, Cells.Imprecise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
set(gca,'YLim', [0 .4]);
set(gca,'XLim', [0 35]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Sago_Validation_TimeCourse.fig');
figure;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, dEdt_Norm);
set(gca,'YLim', [0 .120]);
set(gca,'XLim', [0 35]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Sago_dEdt_Norm.fig');

%% Functions for evaluating indels
function dEdt_Avg = dEdtCalculator(k_total_edit, Time, SimData, Names)
    N_Datapoints = length(Time);
    Total_Cells = SimData(:,1) + SimData(:,2) + SimData(:,3) + SimData(:,4) + SimData(:,5) + SimData(:,6);
    Unedited_Cells = SimData(:,1) + SimData(:,2);
    Norm_dEdtFlux = k_total_edit .* Unedited_Cells .* SimData(:,7) ./ Total_Cells;
    dEdt_Avg = trapz(Time, Norm_dEdtFlux)/max(Time);
end