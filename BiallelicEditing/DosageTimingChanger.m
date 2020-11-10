clear all
%% Script for comparing the effects of dosing with either 6 doses, or a single bolus at a bunch of ages

%% Initialize the model first - Load + Stoptime + Doses
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj');
DoseObj = getdose(m1);
csObj = getconfigset(m1);

TotalDose = sum(DoseObj(2).Amount);
TotalDosemgkg = TotalDose/6.023e23*160e3*1e3/3.5; % (224.5 mg/kg total dose)
SeqDosemgkg = DoseObj(2).Amount(1)/6.023e23*160e3*1e3/3.5; % (23.9 mg/kg)

%% Calculate healing with the doses as set to represent baseline
[T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(2) DoseObj(3)]);
Cbar(1) = C(end);
Sbar(1) = S(end);
Lbar(1) = L(end);

%% Make new dose object
DoseObj(5) = adddose(m1,'s1mplex_1 Bolus','schedule');
DoseObj(6)= adddose(m1,'s1mplex_2 Bolus','schedule');

DoseObj(5).TargetName = 'Liver.S1mplex_1';
DoseObj(6).TargetName = 'Liver.S1mplex_2';

DoseObj(5).Amount = TotalDose;
DoseObj(6).Amount = TotalDose;

DoseObj(5).AmountUnits = 'molecule';
DoseObj(6).AmountUnits = 'molecule';

DoseObj(5).Rate = TotalDose;
DoseObj(6).Rate = TotalDose;

DoseTimes = [0 1 2 3 6 9] * 30;
j = 1;
for i=DoseTimes
    DoseObj(5).Time = DoseTimes(j);
    DoseObj(6).Time = DoseTimes(j);
    [T C S L] = HealingCalculator(m1, csObj(1), [DoseObj(5) DoseObj(6)]);
    Cbar(j+1) = C(end);
    Sbar(j+1) = S(end);
    Lbar(j+1) = L(end);
    j = j+1;
end

bar(Cbar);
set(gca,'YLim',[0 1]);
saveas(gca,'FigureFiles\DoseTiming_C.fig');

figure;
bar(Sbar);
set(gca,'YLim',[0 1]);
saveas(gca,'FigureFiles\DoseTiming_S.fig')

figure;
bar(Lbar);
set(gca,'YLim',[0 1]);
saveas(gca,'FigureFiles\DoseTiming_L.fig')