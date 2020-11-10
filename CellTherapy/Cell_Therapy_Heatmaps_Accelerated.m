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

%% Accelerate the model
CCModel = m1;
sbioaccelerate(CCModel, csObj(1), CCModel_dose);

%% Setup the variables
Granularity = 101;
Progenitor_Dose = 100e6 * linspace(5, 0, Granularity);
Half_Life = linspace(50, 300, Granularity);
GAA_Boost = linspace(0, 1.5, Granularity);

%% Run the files

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose(i), 0.75, Half_Life(j));
    (i-1)*Granularity + j
    "HeatMap 1"
    end
end

Endpoint_Dose_v_Half_Life = Endpoint;

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose(i), GAA_Boost(j), 200);
    (i-1)*Granularity + j
    "HeatMap 2"
    end
end

Endpoint_Dose_v_GAA_Boost = Endpoint;

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, 250e6, GAA_Boost(Granularity + 1 - i), Half_Life(j));
    (i-1)*Granularity + j
    "HeatMap 3"
    end
end

Endpoint_GAA_Boost_v_Half_Life = Endpoint;

save('Cell_Therapy_Heatmap_data_101','Endpoint_Dose_v_GAA_Boost', 'Endpoint_Dose_v_Half_Life', 'Endpoint_GAA_Boost_v_Half_Life', 'Progenitor_Dose', 'GAA_Boost', 'Half_Life');

heatmap(Endpoint_Dose_v_Half_Life(:,:,1));
figure;
heatmap(Endpoint_Dose_v_GAA_Boost(:,:,1));
figure;
heatmap(Endpoint_GAA_Boost_v_Half_Life(:,:,1));
