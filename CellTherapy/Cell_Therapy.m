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

%% Set simulation time to 2 years
csObj = getconfigset(CCModel);
set(csObj(1), 'StopTime', 456);
set(csObj(1), 'TimeUnits','day');

%% Set loss factors 0.02 for heart and 0.008 for muscle
CCModel.Parameters(24).Value = 0.02;
CCModel.Parameters(25).Value = 0.008;

%% Make Progenitor dose object
d_prog = adddose(CCModel, 'Progenitor', 'schedule');
d_prog.TargetName = 'P_Double_Precise';
d_prog.Amount = Progenitor_Dose;
d_prog.AmountUnits = 'cell';
d_prog.Time = 90;
d_prog.TimeUnits = 'day';
d_prog.Rate = d_prog.Amount/24;
d_prog.RateUnits = 'cell/day';

%% Change death rate parameter solely for the added in cells
Cycling_Time = Half_Life;
k_implant_death = 2 * log(2)/Cycling_Time - CCModel.Parameters(26).Value;
Cell_Death_Parameter = addparameter(m1, 'k_implant_death', k_implant_death, 'ValueUnits', '1/day');

Reagent_Name = "P_Double_Precise";
Reaction_Scheme = "P_Double_Precise -> null";
Reaction_Rate = strcat("k_implant_death"," * ",Reagent_Name);
    
Reaction_Scheme_char = convertStringsToChars(Reaction_Scheme);
Reaction_Rate_char = convertStringsToChars(Reaction_Rate);
Reaction_Object = addreaction(m1, Reaction_Scheme_char, 'ReactionRate', Reaction_Rate_char);

Reagent_Name = "M_Double_Precise";
Reaction_Scheme = "M_Double_Precise -> null";
Reaction_Rate = strcat("k_implant_death"," * ",Reagent_Name);
    
Reaction_Scheme_char = convertStringsToChars(Reaction_Scheme);
Reaction_Rate_char = convertStringsToChars(Reaction_Rate);
Reaction_Object = addreaction(m1, Reaction_Scheme_char, 'ReactionRate', Reaction_Rate_char);


%% Retrieve Current dosage information
All_Dose = getdose(CCModel);
CCModel_dose = [All_Dose(5)] ;

%% Accelerate the model
sbioaccelerate(CCModel, csObj(1), CCModel_dose);

