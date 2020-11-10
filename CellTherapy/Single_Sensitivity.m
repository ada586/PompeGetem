%% Thinge to do - 20181212
clear all;
%% Sensitivity Analysis for double edit
% Capture Basal Performance
sbioloadproject('Pompe_Model_v26.sbproj');
CCmodel = m1;
csObj = getconfigset(CCmodel);
set(csObj(1), 'StopTime', 365);
set(csObj(1), 'TimeUnits','day');
All_Dose = getdose(CCmodel);
CCModel_dose = [All_Dose(2) All_Dose(3)];
CCModeldata = sbiosimulate(CCmodel, csObj(1), CCModel_dose);
N_datapoints = size(CCModeldata.Time,1);
BasalHealthyPercent = CCModeldata.Data(N_datapoints, 37)*100;
% All parameters except Healthy ratio are varied 
CyclingParamIndices = [1 2 14 13 11 19 12 7 8 9 10 18]; [7 8 10 11 12 13 14 15 16 18 19];
Perturb = 1e-4;

CCModel_dose = [All_Dose(2) All_Dose(3)];
j = 1;
for i = CyclingParamIndices
    % Load Project
    sbioloadproject('Pompe_Model_v26.sbproj');
    CCmodel = m1;
    csObj = getconfigset(CCmodel);
    set(csObj(1), 'StopTime', 365);
    set(csObj(1), 'TimeUnits','day');
    kvalue = CCmodel.Parameters(i).Value;
    CCmodel.Parameters(i).Value = (1+Perturb) * kvalue;
    CCModeldata = sbiosimulate(CCmodel, csObj(1), CCModel_dose);
    N_datapoints = size(CCModeldata.Time,1);
    PosHealthyPercentSens(j) = CCModeldata.Data(N_datapoints, 37)*100;
    Sensitivity_Axis_Index(j) = convertCharsToStrings(CCmodel.Parameters(i).Name);
    CCmodel.Parameters(i).Value = kvalue;
    j = j+1;
    i
end
Sensitivity_Axis_Index = strrep(Sensitivity_Axis_Index, "_", " ");
Sensitivity = (PosHealthyPercentSens - BasalHealthyPercent)./BasalHealthyPercent/Perturb;

Sens_Bar = bar(Sensitivity);
xtickangle(90)
set(gca, 'fontsize', 30)
xticklabels([" Growth" "Differentiation" "Efficiency" "Precision" "Genome Editor Decay" "Progenitor Affinity" "Cross Correction" "Cellular GAA Degradation" "Serum GAA Decay" "GAA Production" "Double Correction GAA Increase"])

