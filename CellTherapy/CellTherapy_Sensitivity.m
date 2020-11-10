clear all
%% Script to add units of cell and simplex and GAA
cell_unitObj = sbiounit('cell','molecule');
sbioaddtolibrary(cell_unitObj);

s1mplex_unitObj = sbiounit('s1mplex','molecule');
sbioaddtolibrary(s1mplex_unitObj);

GAA_unitObj = sbiounit('GAA','molecule');
sbioaddtolibrary(GAA_unitObj);
%% Script to add dose of precisely edited progenitors or mature cells to the system

%% Load SimBiology Project
sbioloadproject("Pompe_Model_v30_CellDeath.sbproj")
CCModel = m1;

%% Set simulation time to 1 year 3 months
csObj = getconfigset(CCModel);
set(csObj(1), 'StopTime', 456);
set(csObj(1), 'TimeUnits','day');

%% Set loss factors 0.02 for heart and 0.008 for muscle
CCModel.Parameters(24).Value = 0.02;
CCModel.Parameters(25).Value = 0.008;

%% Make Progenitor dose object
d_prog = adddose(CCModel, 'Progenitor', 'schedule');
d_prog.TargetName = 'P_Double_Precise';
d_prog.Amount = 100e6;
d_prog.AmountUnits = 'cell';
d_prog.Time = 90;
d_prog.TimeUnits = 'day';
d_prog.Rate = d_prog.Amount/24;
d_prog.RateUnits = 'cell/day';

%% Change death rate parameter solely for the added in cells
Cycling_Time = 175;
k_implant_death = 2 * log(2)/Cycling_Time - CCModel.Parameters(26).Value;
Cell_Death_Parameter = addparameter(m1, 'k_implant_death', k_implant_death, 'ValueUnits', '1/day');

Reagent_Name = "P_Double_Precise";
Reaction_Scheme = "P_Double_Precise -> null";
Reaction_Rate = strcat("k_implant_death"," * ",Reagent_Name);
    
Reaction_Scheme_char = convertStringsToChars(Reaction_Scheme);
Reaction_Rate_char = convertStringsToChars(Reaction_Rate);
Reaction_Object = addreaction(m1, Reaction_Scheme_char, 'ReactionRate', Reaction_Rate_char);

%% Retrieve Current dosage information
All_Dose = getdose(CCModel);
CCModel_dose = [All_Dose(5)] ;

Progenitor_Dose = 250e6;
GAA_Boost = 0.75;
Half_Life = 175;

% Maturation_Initial = CCModel.Parameters(2).Value
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)
% 
% CCModel.Parameters(2).Value = Maturation_Initial * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% Growth_Initial = CCModel.Parameters(1).Value
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)
% 
% CCModel.Parameters(1).Value = Growth_Initial * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% Progenitor_Dose = Progenitor_Dose * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% Half_Life = Half_Life * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% GAA_Boost = GAA_Boost * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% Cross_Correction = CCModel.Parameters(7).Value
% Cross_Correction = Cross_Correction * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

% Cell_Cross_Correction = CCModel.Parameters(8).Value
% Cell_Cross_Correction = Cell_Cross_Correction * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)

Serum_Cross_Correction = CCModel.Parameters(9).Value
CCModel.Parameters(9).Value = Serum_Cross_Correction * (1+1e-3);
Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
Baseline_Healing(1)

% GAA_Production = CCModel.Parameters(10).Value
% CCModel.Parameters(10).Value = GAA_Production * (1+1e-3);
% Baseline_Healing = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life);
% Baseline_Healing(1)


