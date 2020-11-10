clear all
%% Script for generating both cardiac healing graph, liver healing graph and phenotypes graph with ERT

sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')

% Load the dose and configset
DoseObj = getdose(m1);
csObj = getconfigset(m1);

%% Generating ERT Timeline
% Change loss factors to match ERT
m1.Parameters(24).Value = m1.Parameters(24).Value/4;
m1.Parameters(25).Value = m1.Parameters(25).Value/4;
% Run the model
[T C S L] = HealingCalculator(m1, csObj(1), DoseObj(1));
ERT_Data_Structure.T = T;
ERT_Data_Structure.C = C;
ERT_Data_Structure.S = S;
ERT_Data_Structure.L = L;
C_Fig = plot(T/30,C*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
ERT_C = trapz(T, C)/T(end);
ERT_Data_Structure.ERT_C = ERT_C;
saveas(gca,'FigureFiles\ERT_C_Fig.fig')
figure;
S_Fig = plot(T/30,S*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
ERT_S = trapz(T, S)/T(end);
ERT_Data_Structure.ERT_S = ERT_S;
saveas(gca,'FigureFiles\ERT_S_Fig.fig')
figure;
L_Fig = plot(T/30,L*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
ERT_L = trapz(T, L)/T(end);
ERT_Data_Structure.ERT_L = ERT_L;
saveas(gca,'FigureFiles\ERT_L_Fig.fig')
figure;

%% Reset the loss factors and the dose and generate the same figures for Gene therapy
% Change loss factors to match ERT
m1.Parameters(24).Value = m1.Parameters(24).Value * 4;
m1.Parameters(25).Value = m1.Parameters(25).Value * 4;

[T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
GT_Data_Structure.T = T;
GT_Data_Structure.C = C;
GT_Data_Structure.S = S;
GT_Data_Structure.L = L;
C_Fig = plot(T/30,C*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
GT_C = C(end);
GT_Data_Structure.GT_C = GT_C;
saveas(gca,'FigureFiles\GT_C_Fig.fig')
figure;
S_Fig = plot(T/30,S*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
GT_S = S(end);
GT_Data_Structure.GT_S = GT_S;
saveas(gca,'FigureFiles\GT_S_Fig.fig')
figure;
L_Fig = plot(T/30,L*100);
set(gca,'XLim',[0 12],'YLim',[0 100])
GT_L = L(end);
GT_Data_Structure.GT_L = GT_L;
saveas(gca,'FigureFiles\GT_L_Fig.fig')
figure;

%% Run Gene therapy model and generate the bar chart for the various phenotypes
[T D N] = sbiosimulate(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
Unedited = D(end,1) + D(end,2) + D(end,19) + D(end, 20);
Precise_1_Allele = sum(D(end,3:6))+sum(D(end,15:18));
Imprecise_1_Allele = sum(D(end,7:10))+sum(D(end,21:24));
Precise_2_Allele = sum(D(end,11:12));
Imprecise_2_Allele = sum(D(end,13:14)) + sum(D(end,25:26));
bar([Unedited Precise_1_Allele Imprecise_1_Allele Precise_2_Allele Imprecise_2_Allele])
set(gca,'YScale','log')
set(gca,'YLim', [1e7 1e12]);
saveas(gca,'FigureFiles\GT_Phenotypes.fig')
figure;

