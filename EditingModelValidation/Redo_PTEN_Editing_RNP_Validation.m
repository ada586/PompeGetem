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

[T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
figure;
plot(T, (D(:,5)+D(:,6))/(D(:,1)+D(:,2)+D(:,5)+D(:,6)));

DoseObj(1).Amount = Cas9_mRNA_molecule;
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(1));
figure;
plot(T, (D(:,5)+D(:,6))/(D(:,1)+D(:,2)+D(:,5)+D(:,6)));

