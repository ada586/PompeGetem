function [k_db] = Double_Edit_Rate_Constant(k_Total_Edit,StopTime)
%% Initial concentrations of components to be assigned to specific reactions later
P   = 6.14e10;
S1  = 1;
P1 = 0;
S2  = 1;
P12 = 0;
k1  = k_Total_Edit;
k2  = k_Total_Edit;
ksd1= 0.693;
ksd2= 0.693;
kCCC= 1;

%% Adding the units of cell and s1mplex to the sbioroot
cell_unitObj = sbiounit('cell','molecule');
sbioaddtolibrary(cell_unitObj);

s1mplex_unitObj = sbiounit('s1mplex','molecule');
sbioaddtolibrary(s1mplex_unitObj);

%% Simbiology model for sequential editing
Seq_Edit = sbiomodel("Model for sequential biallelic gene editing");

% Add appropriate compartment for the model
Liver_compartment = addcompartment(Seq_Edit, 'Liver');

% Initialize Species
P_Species_Object = addspecies(Seq_Edit,'P',P, 'InitialAmountUnits', 'cell');
S1_Species_Object = addspecies(Seq_Edit, 'S1', S1, 'InitialAmountUnits', 's1mplex');
P1_Species_Object = addspecies(Seq_Edit,'P1',P1, 'InitialAmountUnits', 'cell');
S2_Species_Object = addspecies(Seq_Edit, 'S2', S2, 'InitialAmountUnits', 's1mplex');
P12_Species_Object = addspecies(Seq_Edit,'P12',P12, 'InitialAmountUnits', 'cell');

% Add parameters
rxn_const_1 = addparameter(Seq_Edit, 'k1', k1, 'ValueUnits', '1/day/s1mplex');
rxn_const_2 = addparameter(Seq_Edit, 'k2', k2, 'ValueUnits', '1/day/s1mplex');
rxn_const_3 = addparameter(Seq_Edit, 'ksd1', ksd1, 'ValueUnits', ' 1/day');
rxn_const_4 = addparameter(Seq_Edit, 'ksd2', ksd2, 'ValueUnits', ' 1/day');

% Reactions
P_2_P1      = addreaction(Seq_Edit, 'P + S1 -> P1');   P_2_P1.ReactionRate = 'k1 * P * S1';
P1_2_P12    = addreaction(Seq_Edit, 'P1 + S2 -> P12');  P1_2_P12.ReactionRate = 'k2 * P1 * S2';
S1_decay    = addreaction(Seq_Edit, 'S1 -> null', 'ReactionRate', 'ksd1 * S1');
S2_decay    = addreaction(Seq_Edit, 'S2 -> null', 'ReactionRate', 'ksd2 * S2');

%% Simbiology Model for simultanous editing using k.S harmonic means
% Define the harmonic mean rate constant
k12 = kCCC * k1*k2/(k1+k2)*1/(S1+S2);
% Setup Model
SimulEdit_HM = sbiomodel('Model for simultaneous biallelic editing using k.S Harmonic Means');

% Setup Compartment
Liver_compartment = addcompartment(SimulEdit_HM, 'Liver');

% Initialize Species
P_Species_Object = addspecies(SimulEdit_HM, 'P', P, 'InitialAmountUnits', 'cell');
S1_Species_Object = addspecies(SimulEdit_HM, 'S1', S1, 'InitialAmountUnits', 's1mplex');
S2_Species_Object = addspecies(SimulEdit_HM, 'S2', S2, 'InitialAmountUnits', 's1mplex');
P12_Species_Object = addspecies(SimulEdit_HM, 'P12', P12, 'InitialAmountUnits', 'cell');

% Add parameters
rxn_const_1 = addparameter(SimulEdit_HM, 'k12', k12, 'ValueUnits', '1/day/s1mplex/s1mplex');
rxn_const_3 = addparameter(SimulEdit_HM, 'ksd1', ksd1, 'ValueUnits', ' 1/day');
rxn_const_4 = addparameter(SimulEdit_HM, 'ksd2', ksd2, 'ValueUnits', ' 1/day');

% Add reactions
P_2_P12     = addreaction(SimulEdit_HM, 'P + S1 + S2 -> P12', 'ReactionRate', 'k12 * S1 * S2 * P');
S1_decay    = addreaction(SimulEdit_HM, 'S1 -> null', 'ReactionRate', 'ksd1 * S1');
S2_decay    = addreaction(SimulEdit_HM, 'S2 -> null', 'ReactionRate', 'ksd2 * S2');

%% Set solving conditions and calculation of residuals
% Setting a stop time of 7 days for simulation of Seq_Edit
% Capturing current simulation settings
New_k12 = SimulEdit_HM.Parameters(1).Value;

csObj = addconfigset(Seq_Edit, 'Active');
csObj.StopTime = StopTime;
csObj.TimeUnits = 'day';

HM_csObj = addconfigset(SimulEdit_HM, 'Active');
HM_csObj.StopTime = StopTime;
HM_Obj.TimeUnits = 'day';

for i=1:5
    % Run the simulation
    [Time SimData Names] = sbiosimulate(Seq_Edit, csObj);
    
    % Extract the P12 and P data from the SimData
    P_Data = SimData(:,1);
    P12_Data = SimData(:,5);
    
    % Set conditions for solving with HM_1
    SimulEdit_HM.Parameters(1).Value = New_k12;
    
    [HM_Time SimData Names] = sbiosimulate(SimulEdit_HM, HM_csObj);
    
    % Extract Data for HM_1
    P_HM = SimData(:,1);
    P12_HM = SimData(:,4);
    
    %plot(Time, P12_Data, Time, P12_HM);
    % Extract old k_12 value, modify if by multiplying P_12/P_HM
    New_k12 = New_k12 * P12_Data(length(Time))/P12_HM(length(HM_Time));
    i;
    P12_Data(length(Time))/P12_HM(length(HM_Time));
end
    k_db = New_k12;
end

