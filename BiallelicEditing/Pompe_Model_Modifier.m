clear all
%% Script to generate an updated Pompe Model sbproj file
sbioloadproject('../EditingModelValidation/SimBiologyProjectFiles/PTEN_Editing_Project_Wei_2020.sbproj')
ms_RNP_Total = m1.Species(7).InitialAmount;
ms_Cells_Total = m1.Species(1).InitialAmount + m1.Species(2).InitialAmount;
RNP_per_Cell = ms_RNP_Total/ms_Cells_Total;

sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v27.sbproj');
hu_Cells_Total = m1.Species(1).InitialAmount + m1.Species(2).InitialAmount;
hu_RNP_Total = RNP_per_Cell * hu_Cells_Total;

load('DataMatlabFiles\EditingRateConstants.mat')

csObj = getconfigset(m1);
DoseObj = getdose(m1);

%% Simulate ERT using relevant loss factors
m1.Parameters(25).Value = 0.008 / 4; % For ERT set to 0.002, 0.008 normal
m1.Parameters(24).Value = 0.02 / 4;  % For ERT set to 0.005, 0.02 normal

[T C S L] = HealingCalculator(m1, csObj(1), DoseObj(1));
Desired_C = trapz(T,C)/max(T);
Desired_S = trapz(T,S)/max(T);
Desired_L = trapz(T,L)/max(T);
plot(T,C);

%% Iterate to find the acceptable dose
% Reset the parameter values for loss factors
m1.Parameters(25).Value = 0.008; % For ERT set to 0.002, 0.008 normal
m1.Parameters(24).Value = 0.02;  % For ERT set to 0.005, 0.02 normal
m1.Parameters(15).Value = k_Range(3);
m1.Parameters(16).Value = BiAllelic_k_Range(3);

%% Start the iteration
DoseAmount_Old = DoseObj(2).Amount;
m1.Parameters(19).Value = 0.5;

for i=1:20
    [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
    Healing_C = C(end);
    Modifier(i) = Desired_C/Healing_C;
    DoseObj(2).Amount = DoseObj(2).Amount * Modifier(i);
    DoseObj(3).Amount = DoseObj(3).Amount * Modifier(i);
end

DoseAmount_New = DoseObj(2).Amount;

Dose_per_kg = DoseAmount_New(1) / 6.023e23 * 160e3 * 1e3/3.5 ; %mg/kg for newborn male infant

csObj(1).StopTime = 30;
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(2));
m1.Parameters(15).Value = k_Range(3);
m1.Parameters(16).Value = BiAllelic_k_Range(3);
EditRate = m1.Parameters(15).Value * (D(:,1) + D(:,2)) .* D(:,27);
TotalEdits = trapz(T,EditRate);
TotalEditPercentage = TotalEdits / sum(D(end,1:26))*100;

%sbiosaveproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj','m1')