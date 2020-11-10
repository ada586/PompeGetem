clear all
%% Script to add units of cell and simplex and GAA
cell_unitObj = sbiounit('cell','molecule');
sbioaddtolibrary(cell_unitObj);

s1mplex_unitObj = sbiounit('s1mplex','molecule');
sbioaddtolibrary(s1mplex_unitObj);

GAA_unitObj = sbiounit('GAA','molecule');
sbioaddtolibrary(GAA_unitObj);

%% Load Simbiology Project
sbioloadproject("Pompe_Model_v26.sbproj")
CCModel = m1;

%% Set simulation time to 1 year
csObj = getconfigset(CCModel);
set(csObj(1), 'StopTime', 365);
set(csObj(1), 'TimeUnits','day')
All_Dose = getdose(CCModel);

%% Setting up doses for the simulation
% All_Dose(1) = biweekly rhGAA, ERT
% All_Dose(2) = 6 x 614 mg/kg S1mplex 1
% All_Dose(3) = 6 x 614 mg/kg S1mplex 2
% All_Dose(4) = 6 x 1228 mg/kg S1mplex 1
% Defining a dose amount using linspace
% ModelDose_Amount = 1.1407e17*linspace(3.4, 7.233, 6); %30.3 mg/kg
% ModelDose_Rate = 1.1407e17*linspace(3.4, 7.233, 6); %30.3 mg/kg

ModelDose_Amount = 9.7899e16*linspace(3.4, 7.233, 6); %30.3 mg/kg
ModelDose_Rate = 9.7899e16*linspace(3.4, 7.233, 6); %30.3 mg/kg

All_Dose(2).Amount = ModelDose_Amount;
All_Dose(3).Amount = ModelDose_Amount;
All_Dose(4).Amount = ModelDose_Amount * 2;
All_Dose(2).Rate = ModelDose_Amount;
All_Dose(3).Rate = ModelDose_Amount;
All_Dose(4).Rate = ModelDose_Amount * 2;
CCModel_dose = [All_Dose(2) All_Dose(3)];

%% Simulate Model
CCModeldata = sbiosimulate(CCModel, csObj(1), CCModel_dose);
Model_Dose_perkg = 1.1404;

for i=1:10
    ModelDose_Amount = Model_Dose_perkg*linspace(3.4, 7.233, 6); %30.3 mg/kg
    ModelDose_Rate = Model_Dose_perkg*linspace(3.4, 7.233, 6); %30.3 mg/kg
    
    All_Dose(2).Amount = ModelDose_Amount;
    All_Dose(3).Amount = ModelDose_Amount;
    All_Dose(4).Amount = ModelDose_Amount * 2;
    All_Dose(2).Rate = ModelDose_Amount;
    All_Dose(3).Rate = ModelDose_Amount;
    All_Dose(4).Rate = ModelDose_Amount * 2;
    CCModel_dose = [All_Dose(2) All_Dose(3)];
    
    CCModeldata = sbiosimulate(CCModel, csObj(1), CCModel_dose);
    Time = CCModeldata.Time;
    Total_C = CCModeldata.Data(:,33) + CCModeldata.Data(:,34);
    Fraction_C = CCModeldata.Data(:,34)./Total_C * 100;
    Healing = Fraction_C(length(Time))
    i
    Model_Dose_perkg = Model_Dose_perkg/Healing*32.6748;
end


%% Calculate percentage of normal phenotype cells
Time = CCModeldata.Time;

% Liver Healing
NormalPercent = CCModeldata.Data(:,37)*100;
plot(CCModeldata.Time,CCModeldata.Data(:,37)*100);
ylim([0 100]);
title("Liver")
figure;
trapz(CCModeldata.Time,CCModeldata.Data(:,37)*100)/365;

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
