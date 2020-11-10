clear all
%% Script to generate the distribution of phenotypes bar chart
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')
DoseObj = getdose(m1);
csObj = getconfigset(m1);

%% Run Gene therapy model and generate the bar chart for the various phenotypes
[T D N] = sbiosimulate(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
Unedited = D(end,1) + D(end,2) + D(end,19) + D(end, 20);
Precise_1_Allele = sum(D(end,3:6));
Imprecise_1_Allele = sum(D(end,7:10))+sum(D(end,21:24));
Precise_2_Allele = sum(D(end,11:12));
Precise_1_Imprecise_1 = sum(D(end,15:18));
Imprecise_2_Allele = sum(D(end,13:14)) + sum(D(end,25:26));
PhenotypesBarData = [Unedited Precise_1_Allele Imprecise_1_Allele Precise_2_Allele Precise_1_Imprecise_1 Imprecise_2_Allele];
bar(PhenotypesBarData)
set(gca,'YScale','log')
set(gca,'YLim', [1e7 1e12]);
saveas(gca,'FigureFiles\GT_Phenotypes.fig');

%% Graphing Distance calculator
DecadeHeight = (354.1936 - 293.7434)/(12 - 7);

BarHeights = (log10(PhenotypesBarData) - 7) * DecadeHeight;