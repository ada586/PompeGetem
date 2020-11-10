clear all

%% Sensitivity analysis of the parameters

sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')

% Load the dose and configset
DoseObj = getdose(m1);
csObj = getconfigset(m1);
ParamList = [1 2 14 11 9 10 20 21 7 8 13 19 18 12];
j = 1;
for i = ParamList
    SensitivityParameter(j).Name = m1.Parameters(i).Name;
    j = j+1;
end

csObj(1).StopTime = 365;
j = 1;
Perturb = 1e-4;

[T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
C_old = C(end);
S_old = S(end);
L_old = L(end);
    
for i = ParamList
    k_value = m1.Parameters(i).Value;
    m1.Parameters(i).Value = (1+Perturb) * k_value;
    [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
    C_new = C(end);
    S_new = S(end);
    L_new = L(end);
    SensitivityParameter(j).C_Sensitivity = (C_new - C_old)/C_old/Perturb;
    SensitivityParameter(j).S_Sensitivity = (S_new - S_old)/S_old/Perturb;
    SensitivityParameter(j).L_Sensitivity = (L_new - L_old)/L_old/Perturb;
    m1.Parameters(i).Value = k_value;
    Sensitivity_Axis_Title(j) = convertCharsToStrings(m1.Parameters(i).Name);
    Sensitivity_Axis_Title(j) = strrep(Sensitivity_Axis_Title(j),"_"," ");
    Sensitivity_Value(j) = SensitivityParameter(j).C_Sensitivity;
    j=j+1;
end
% 
barh(fliplr(abs(Sensitivity_Value(1:2))))
yticklabels(fliplr(Sensitivity_Axis_Title(1:2)))
%xtickangle(45)
saveas(gca, 'SensitivityFigures\v_27_Morphogenesis_Perturb_4.fig');
figure;

barh(fliplr(abs(Sensitivity_Value(3:end))))
yticklabels(fliplr(Sensitivity_Axis_Title(3:end)))
%xtickangle(45)
saveas(gca, 'SensitivityFigures\v_27_Others_Perturb_4.fig');