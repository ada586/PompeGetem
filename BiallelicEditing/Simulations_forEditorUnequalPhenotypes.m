clear all
%% Script to run experiments required for reviewer 4

%% Load project, configset and doses
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')
csObj = getconfigset(m1);
DoseObj = getdose(m1);
%% Adding in Precise Editing Modifier to make Allele 1 editing less precise