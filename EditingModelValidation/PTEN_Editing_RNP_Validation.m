clear all
%% Validation of GETEM for PTEN editing using SpyCas9 RNP
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

% Load the PTEN editing with RNP
sbioloadproject('SimBiologyProjectFiles\PTEN_Editing_Project_Wei_2020.sbproj')
Ms_RNP_Dose = m1.Species(7).InitialAmount;

% k_total_edit = dEdt/RNP, calculate mean lower and upper bounds for
% k_total_edit

k_total_edit = [(Aggregate_dEdt_Mean - Aggregate_dEdt_Std)/Ms_RNP_Dose, ...
    (Aggregate_dEdt_Mean)/Ms_RNP_Dose, ...
    (Aggregate_dEdt_Mean + Aggregate_dEdt_Std)/Ms_RNP_Dose];

for i=1:3
    m1.Parameters(9).Value = k_total_edit(i);
    [Time, SimData, Names] = sbiosimulate(m1);
    N_Datapoints = length(Time);
    TotalCells(i) = SimData(N_Datapoints,1) + SimData(N_Datapoints,2) + SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditedCells(i) = SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditFraction(i) = EditedCells(i)/TotalCells(i);
end

save('ValidationdEdtData\PTEN_Editing_Predictions.mat','EditFraction');

