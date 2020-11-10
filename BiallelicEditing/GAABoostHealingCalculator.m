clear all
%% Excess GAA Production over single allele calculator

%% Initialize the model and the dose objects
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')
DoseObj = getdose(m1);
csObj = getconfigset(m1);

%% Increase the double allele GAA boost and check
GAABoostArray = [0 0.25 0.5 0.75 1 1.25 1.50];

for i=1:7
    m1.Parameters(18).Value = GAABoostArray(i);
    [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
    Cbar(i) = C(end);
    Sbar(i) = S(end);
    Lbar(i) = L(end);
end
