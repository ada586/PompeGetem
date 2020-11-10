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
csObj(1).StopTime = 30;
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(2));
EditRate = m1.Parameters(15).Value * (D(:,1) + D(:,2)) .* D(:,27);
TotalEdits = trapz(T,EditRate);
TotalEditPercentage = TotalEdits / sum(D(end,1:26))*100;

% Adjusting total editing rate
m1.Parameters(15).Value = m1.Parameters(15).Value/TotalEditPercentage;
m1.Parameters(16).Value = Double_Edit_Rate_Constant(m1.Parameters(15).Value, 35);

% Reset Stoptime
csObj(1).StopTime = 365;
% Accelerate the model
sbioaccelerate(m1, csObj(1), [DoseObj(2) DoseObj(3)]);



% Generate Efficiency vs. Precision
Efficiency_linspaced = linspace(20, 0, Granularity);
%Efficiency_linspaced = logspace(log10(20), log10(0.5), Granularity);
P_IP_Percent  = linspace(10, 90, Granularity);
P_IP_linspaced = P_IP_Percent ./ (100 - P_IP_Percent);

for i=1:Granularity
    m1.Parameters(14).Value = Efficiency_linspaced(i);
    for j=1:Granularity
        m1.Parameters(13).Value = P_IP_linspaced(j);
        [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
        C_Heatmap(i,j) = C(end);
        S_Heatmap(i,j) = S(end);
        L_Heatmap(i,j) = L(end);
        [i j]
    end
end

C_h = heatmap(C_Heatmap);
% C_h.YDisplayLabels = Efficiency_linspaced;
% C_h.XDisplayLabels = P_IP_Percent;
C_h.ColorLimits = [0 1];
C_h.Colormap = Colormaps.C.Array;
C_h.GridVisible = 'off';
C_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\Efficiency_Precision_C.fig')
figure;

S_h = heatmap(S_Heatmap);
% S_h.YDisplayLabels = Efficiency_linspaced;
% S_h.XDisplayLabels = P_IP_Percent;
S_h.ColorLimits = [0 1];
S_h.Colormap = Colormaps.S.Array;
S_h.GridVisible = 'off';
S_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\Efficiency_Precision_S.fig')
figure;

L_h = heatmap(L_Heatmap);
% S_h.YDisplayLabels = Efficiency_linspaced;
% S_h.XDisplayLabels = P_IP_Percent;
L_h.ColorLimits = [0 1];
L_h.Colormap = Colormaps.L.Array;
L_h.GridVisible = 'off';
L_h.CellLabelColor = 'none';
saveas(gca, 'FigureFiles\Efficiency_Precision_L.fig')

save('DataMatlabFiles\EfficiencyPrecisionHeatmaps.mat','C_Heatmap', 'S_Heatmap', 'L_Heatmap', 'Efficiency_linspaced', 'P_IP_linspaced');
