clear all
%% Calculation of the range of total editing rate constants and corresponding percentage editing numbers from the validation dataset

%% Calculation of the ms RNP per cell from the PTEN validation Simbiology model
sbioloadproject('../EditingModelValidation/SimBiologyProjectFiles/PTEN_Editing_Project_Wei_2020.sbproj')
ms_RNP_Total = m1.Species(7).InitialAmount;
ms_Cells_Total = m1.Species(1).InitialAmount + m1.Species(2).InitialAmount;
RNP_per_Cell = ms_RNP_Total/ms_Cells_Total;

%% Calculation of range of dE_dt from validation data
load('../EditingModelValidation/ValidationdEdtData/ValidationStructureFile_dEdt_Norm.mat');
for i=1:4
    dEdt_Data(i).Mean = mean(ValidationModel(i).NormalizeddEdt);
    dEdt_Data(i).Std  = std(ValidationModel(i).NormalizeddEdt);
end

Aggregate_dEdt_Mean = mean([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
Aggregate_dEdt_Std = std([dEdt_Data(1).Mean, dEdt_Data(2).Mean, dEdt_Data(3).Mean, dEdt_Data(4).Mean]);
dEdt_Range = [Aggregate_dEdt_Mean-Aggregate_dEdt_Std, Aggregate_dEdt_Mean, Aggregate_dEdt_Mean+Aggregate_dEdt_Std];

%% Calculation of the human rate constant from the Pompe model basis
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v26.sbproj');
hu_Cells_Total = m1.Species(1).InitialAmount + m1.Species(2).InitialAmount;
hu_RNP_Total = RNP_per_Cell * hu_Cells_Total;
k_Range = dEdt_Range./hu_RNP_Total;

for i=1:3
    BiAllelic_k_Range(i) = Double_Edit_Rate_Constant(k_Range(i), 35);
end
% Initiate the rest of the Pompe model
csObj = getconfigset(m1);
DoseObj = getdose(m1);

% To fugure out the efficiency of editing from the rate constant and the
% dose calculated above
for i=1:3
    m1.Parameters(15).Value = k_Range(i);
    m1.Parameters(16).Value = BiAllelic_k_Range(i);
    m1.Species(27).InitialAmount = hu_RNP_Total;
    csObj(1).StopTime = 7;
    [Time SimData Names] = sbiosimulate(m1, csObj(1));
    dEdt_Flux = m1.Parameters(15).Value * (SimData(:,1) + SimData(:,2)) .* SimData(:,27);
    Total_Cells = SimData(:,1) + SimData(:,2) + SimData(:,7) + SimData(:,8) + SimData(:,3) + SimData(:,4);
    dEdt_Flux = m1.Parameters(15).Value * (SimData(:,1) + SimData(:,2)) .* SimData(:,27) ./ Total_Cells;
    EditPercent(i) = trapz(Time, dEdt_Flux);
end

% Rescale k to be responsive to 1% Editing
k_Range_Rescaled = k_Range / EditPercent(2) / 100;
for i=1:3
    BiAllelic_k_Range_Rescaled(i) = Double_Edit_Rate_Constant(k_Range(i), 35);
end

save('DataMatlabFiles\EditingRateConstants.mat','k_Range_Rescaled','BiAllelic_k_Range_Rescaled','k_Range','BiAllelic_k_Range');

