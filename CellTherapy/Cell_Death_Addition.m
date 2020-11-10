clear all
% Script for adding in the death rate parameter, adjusting the
% differentiation rate (via a rule) and adding in all of the death
% reactions for even species 2:26.

%%
% Load project with cell deaths built in and initialize a 100 day liver
% turnover death rate
sbioloadproject('Pompe_Model_v27_CellDeath.sbproj');

k_mature_death = log(2)/50;  % 1/day

%% Add parameter for cell death
Cell_Death_Parameter = addparameter(m1, 'k_mature_death', k_mature_death, 'ValueUnits', '1/day');

%% Add rule modifying the death rate
Cell_Death_Rule = addrule(m1, 'm1.Parameters(2).Value = 0.99 * (m1.Parameters(1).Value + m1.Parameters(26).Value)');

%% Add reactions for species 2:26, even numbers only
for i = 1:13
    j = 2 * i;
    Reagent_Name = convertCharsToStrings(m1.Species(j).Name);
    Reaction_Scheme = strcat(Reagent_Name, " -> null");
    Reaction_Rate = strcat(convertCharsToStrings(m1.Parameters(26).Name)," * ",Reagent_Name);
    
    Reaction_Scheme_char = convertStringsToChars(Reaction_Scheme);
    Reaction_Rate_char = convertStringsToChars(Reaction_Rate);
    
    Reaction_Object(i) = addreaction(m1, Reaction_Scheme_char, 'ReactionRate', Reaction_Rate_char);
end


