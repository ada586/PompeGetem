clear all
%% Script to add units of cell and simplex and GAA
cell_unitObj = sbiounit('cell','molecule');
sbioaddtolibrary(cell_unitObj);

s1mplex_unitObj = sbiounit('s1mplex','molecule');
sbioaddtolibrary(s1mplex_unitObj);

GAA_unitObj = sbiounit('GAA','molecule');
sbioaddtolibrary(GAA_unitObj);

%% Load Simbiology Project
sbioloadproject("Pompe_Model_v29_CellDeath.sbproj")
CCModel = m1;

%% Set death rate;
CCModel.Parameters(26).Value = log(2)/25;

% Equations for altering growth rates to compensate for death rates
% CCModel.Parameters(2).Value = 100 * (CCModel.Parameters(26).Value + log(2.8964)/365);
% CCModel.Parameters(1).Value = 1.01 * CCModel.Parameters(2).Value - CCModel.Parameters(26).Value;

%% Set simulation time to 1 year
csObj = getconfigset(CCModel);
set(csObj(1), 'StopTime', 365);
set(csObj(1), 'TimeUnits','day')
All_Dose = getdose(CCModel);

%% Set Cardiac loss parameter and skeletal loss parameter
% CCModel.Parameters(25).Value = 0.008 / 4; % For ERT set to 0.002, 0.008 normal
% CCModel.Parameters(24).Value = 0.02 / 4;  % For ERT set to 0.005, 0.02 normal
% For ERT use Dose(1), for single s1mplex big dose use Dose(4)
CCModel_dose = [All_Dose(2) All_Dose(3)];

%% Simulate Model
CCModeldata = sbiosimulate(CCModel, csObj(1), CCModel_dose);

%% Calculate percentage of normal phenotype cells
Time = CCModeldata.Time;

% Liver Healing
NormalPercent = CCModeldata.Data(:,37)*100;
plot(CCModeldata.Time,CCModeldata.Data(:,37)*100);
ylim([0 100]);
title("Liver")
figure;
trapz(CCModeldata.Time,CCModeldata.Data(:,37)*100)/365;

%% Calculate Genome Percentages
Last_Row = CCModeldata.Data(length(Time),:);
Unedited_Total = Last_Row(1) + Last_Row(2) + Last_Row(19) + Last_Row(20);
Single_Precise_Edit = Last_Row(3) + Last_Row(4) + Last_Row(5) + Last_Row(6);
Single_Imprecise_Edit = Last_Row(7) + Last_Row(8) + Last_Row(9) + Last_Row(10) + Last_Row(21) + Last_Row(22) + Last_Row(23) + Last_Row(24);
Double_Precise_Edit = Last_Row(11) + Last_Row(12);
Double_Imprecise_Edit = Last_Row(13) + Last_Row(14) + Last_Row(25) + Last_Row(26);
Single_Precise_Single_Imprecise = Last_Row(15) + Last_Row(16) + Last_Row(17) + Last_Row(18);

%% Calculate Efficiency at first month
Last_Row = CCModeldata.Data(191,:);
Unedited_Total_191 = Last_Row(1) + Last_Row(2) + Last_Row(19) + Last_Row(20);
Single_Precise_Edit_191 = Last_Row(3) + Last_Row(4) + Last_Row(5) + Last_Row(6);
Single_Imprecise_Edit_191 = Last_Row(7) + Last_Row(8) + Last_Row(9) + Last_Row(10) + Last_Row(21) + Last_Row(22) + Last_Row(23) + Last_Row(24);
Double_Precise_Edit_191 = Last_Row(11) + Last_Row(12);
Double_Imprecise_Edit_191 = Last_Row(13) + Last_Row(14) + Last_Row(25) + Last_Row(26);
Single_Precise_Single_Imprecise_191 = Last_Row(15) + Last_Row(16) + Last_Row(17) + Last_Row(18);

Efficiency_In_Situ = (Single_Precise_Edit_191+Single_Imprecise_Edit_191)/Unedited_Total_191; %About 2.26% Efficiency 
Dbl_Efficiency_In_Situ = (Double_Precise_Edit_191 + Double_Imprecise_Edit_191 + Single_Precise_Single_Imprecise_191)/Unedited_Total_191; % About 0.0176% Efficient
%% Other plots
% Cardiac Healing
Total_C = CCModeldata.Data(:,33) + CCModeldata.Data(:,34);
Fraction_C = CCModeldata.Data(:,34)./Total_C * 100;
plot(CCModeldata.Time,Fraction_C);
ylim([0 100]);
title("Heart")
figure;
trapz(CCModeldata.Time,Fraction_C)/365;
Fraction_C(length(Time))

% Skeletal Healing
Total_S = CCModeldata.Data(:,35) +CCModeldata.Data(:,36);
Fraction_S = CCModeldata.Data(:,36)./Total_S * 100;
plot(CCModeldata.Time,Fraction_S);
ylim([0 100]);
title("Muscle")
trapz(CCModeldata.Time,Fraction_S)/365;
figure;

% Total Liver Growth
plot(Time,CCModeldata.Data(:,30)+ CCModeldata.Data(:,31))
Total_Liver = CCModeldata.Data(:,30)+ CCModeldata.Data(:,31);