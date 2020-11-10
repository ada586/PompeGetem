clear all
%% Set the granularity
Granularity = 100;
%% Heatmap Colormap maker
Colormaps = ColorMapMaker(Granularity);
%% Script for generating the heatmaps
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')
csObj = getconfigset(m1);
DoseObj = getdose(m1);
% Current total edit rate represents a certain editing efficiency,
% adjusting it to reflect 1% editing
% Run model for 28 days, check percentage of liver edited using editing
% equation flux with a single editor
% csObj(1).StopTime = 30;
% [T D N] = sbiosimulate(m1, csObj(1), DoseObj(2));
% EditRate = m1.Parameters(15).Value * (D(:,1) + D(:,2)) .* D(:,27);
% TotalEdits = trapz(T,EditRate);
% TotalEditPercentage = TotalEdits / sum(D(end,1:26))*100;
% 
% % Adjusting total editing rate
% m1.Parameters(15).Value = m1.Parameters(15).Value/TotalEditPercentage;
% m1.Parameters(16).Value = Double_Edit_Rate_Constant(m1.Parameters(15).Value, 35);

%% Set model up for running heatmap
% Reset Stoptime
csObj(1).StopTime = 365;
% Accelerate the model
sbioaccelerate(m1, csObj(1), [DoseObj(2) DoseObj(3)]);

%% Generate dose vs progenitor affinity heatmap
% Generate dose vector from 1 to 100 mg/kg
DoseLinspaced = logspace(log10(100), log10(10), Granularity);
% DoseLinspaced = linspace(50, 5, Granularity);
DoseLinspacedMolecules = DoseLinspaced * 1e-3 / 160e3 * 6.023e23 * 3.5;
ScalingFactorDose = DoseLinspacedMolecules ./ DoseObj(2).Amount(1);

% ProgenitorAffinityLinspaced = logspace(log10(0.01), log10(1), Granularity);
ProgenitorAffinityLinspaced = linspace(0, 1, Granularity);

OrigDose2Vector = DoseObj(2).Amount;
OrigDose3Vector = DoseObj(3).Amount;

for i=1:Granularity
    DoseObj(2).Amount = ScalingFactorDose(i) * OrigDose2Vector;
    DoseObj(3).Amount = ScalingFactorDose(i) * OrigDose3Vector;
    for j=1:Granularity
        m1.Parameters(19).Value = ProgenitorAffinityLinspaced(j);
        [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
        C_Heatmap(i,j) = C(end);
        S_Heatmap(i,j) = S(end);
        L_Heatmap(i,j) = L(end);
        [i j]
    end
end

C_h = heatmap(C_Heatmap);
C_h.YDisplayLabels = DoseLinspaced;
C_h.XDisplayLabels = ProgenitorAffinityLinspaced;
C_h.ColorLimits = [0 1];
C_h.Colormap = Colormaps.C.Array;
C_h.GridVisible = 'off';
C_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\DoseProgenitorAffinityHeatmaps_C.fig')
figure;

S_h = heatmap(S_Heatmap);
S_h.YDisplayLabels = DoseLinspaced;
S_h.XDisplayLabels = ProgenitorAffinityLinspaced;
S_h.ColorLimits = [0 1];
S_h.Colormap = Colormaps.S.Array;
S_h.GridVisible = 'off';
S_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\DoseProgenitorAffinityHeatmaps_S.fig')
figure;

L_h = heatmap(L_Heatmap);
L_h.YDisplayLabels =DoseLinspaced;
L_h.XDisplayLabels = ProgenitorAffinityLinspaced;
L_h.ColorLimits = [0 1];
L_h.Colormap = Colormaps.L.Array;
L_h.GridVisible = 'off';
L_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\DoseProgenitorAffinityHeatmaps_L.fig')

save('DataMatlabFiles\DoseProgenitorAffinityHeatmaps.mat','C_Heatmap', 'S_Heatmap', 'L_Heatmap', 'DoseLinspaced', 'ProgenitorAffinityLinspaced');
