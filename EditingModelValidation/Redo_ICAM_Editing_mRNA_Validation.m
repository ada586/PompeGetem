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
DoseObj(1).Amount = Cas9_mRNA_molecule;

% Set stop time to 7 days
csObj = getconfigset(m1);
csObj(1).StopTime = 7;

% Set precise editing rate constant to zero
m1.Parameters(8).Value = 0;

%% Sanity check time
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
figure;
plot(T, (D(:,5)+D(:,6))/(D(:,1)+D(:,2)+D(:,5)+D(:,6)));

DoseObj(1).Amount = Cas9_mRNA_molecule;
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
figure;
plot(T, (D(:,5)+D(:,6))/(D(:,1)+D(:,2)+D(:,5)+D(:,6)));

%% Import k values from the mRNA simulation
load('ValidationdEdtData\ValidationStructureFile_dEdt_Norm.mat')
% Use the dEdt data

for i=1:4
    dEdt_Data(i).Mean = mean(ValidationModel(i).NormalizeddEdt);
    dEdt_Data(i).Std  = std(ValidationModel(i).NormalizeddEdt);
end

Aggregate_dEdt_Mean = mean([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
Aggregate_dEdt_Std = std([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);

dEdt_Range = [Aggregate_dEdt_Mean-Aggregate_dEdt_Std Aggregate_dEdt_Mean Aggregate_dEdt_Mean+Aggregate_dEdt_Std];

for i=1:3
    Desired_dEdt = dEdt_Range(i);
    for j=1:100
        k_imprecise_edit = m1.Parameters(9).Value;
        [Time SimData Names] = sbiosimulate(m1, csObj(1), DoseObj(1));
        dEdt_Avg = dEdtCalculator(k_imprecise_edit, Time, SimData, Names);
        k_modifier = Desired_dEdt/dEdt_Avg;
        k_imprecise_edit = k_imprecise_edit * k_modifier;
        [i, j];
    end
    k_total_edit(i) = k_imprecise_edit;
end

for i=1:3
    m1.Parameters(9).Value = k_total_edit(i);
    [T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
    EditedAmount = (D(:,5)+D(:,6))./(D(:,1)+D(:,2)+D(:,5)+D(:,6));
    EditRangeFullDose(i) = EditedAmount(end);
end
% %% Medium Dose - Based on limiting reagent
% ms_wt = (17.8 + 28.9)/2; % g
% Cas9_mRNA_Dose = 1/4 * ms_wt/1000 * 1e-3; %g
% Cas9_mRNA_MW = (1368 * 320.5) + 159; % g/mol
% Cas9_mRNA_mol = Cas9_mRNA_Dose / Cas9_mRNA_MW;
% Cas9_mRNA_molecule = Cas9_mRNA_mol * 6.023e23;
% 
% DoseObj(1).Amount = Cas9_mRNA_molecule;
% for i=1:3
%     m1.Parameters(9).Value = k_total_edit(i);
%     [T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
%     EditedAmount = (D(:,5)+D(:,6))./(D(:,1)+D(:,2)+D(:,5)+D(:,6));
%     EditRangeMedDose(i) = EditedAmount(end);
% end
% %% Low Dose - based on limiting reagent
% ms_wt = (17.8 + 28.9)/2; % g
% Cas9_mRNA_Dose = 1/6 * ms_wt/1000 * 1e-3; %g
% Cas9_mRNA_MW = (1368 * 320.5) + 159; % g/mol
% Cas9_mRNA_mol = Cas9_mRNA_Dose / Cas9_mRNA_MW;
% Cas9_mRNA_molecule = Cas9_mRNA_mol * 6.023e23;
% 
% DoseObj(1).Amount = Cas9_mRNA_molecule;
% for i=1:3
%     m1.Parameters(9).Value = k_total_edit(i);
%     [T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
%     EditedAmount = (D(:,5)+D(:,6))./(D(:,1)+D(:,2)+D(:,5)+D(:,6));
%     EditRangeMedDose(i) = EditedAmount(end);
% end
% DoseObj(1).Amount = Cas9_mRNA_molecule/6;
% for i=1:3
%     m1.Parameters(9).Value = k_total_edit(i);
%     [T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
%     EditedAmount = (D(:,5)+D(:,6))./(D(:,1)+D(:,2)+D(:,5)+D(:,6));
%     EditRangeLowDose(i) = EditedAmount(end);
% end

%% Functions for evaluating indels
function dEdt_Avg = dEdtCalculator(k_total_edit, Time, SimData, Names)
    N_Datapoints = length(Time);
    Total_Cells = SimData(:,1) + SimData(:,2) + SimData(:,3) + SimData(:,4) + SimData(:,5) + SimData(:,6);
    Unedited_Cells = SimData(:,1) + SimData(:,2);
    Norm_dEdtFlux = k_total_edit .* Unedited_Cells .* SimData(:,7) ./ Total_Cells;
    dEdt_Avg = trapz(Time, Norm_dEdtFlux)/max(Time);
end