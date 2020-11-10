%% Function to generate RNP dynamics, Edited Cell Flux Dynamics using an accelerated model
function [Time, RNP, dE_dt_Norm] = RNP_Flux_Dynamics(m1, k_total_edit)
    % Start by initiating the model as written
    csObj = getconfigset(m1, 'Active');
    csObj.TimeUnits = 'day';
    doseObj = m1.getdose;
    
    % Find precision relationship and assign editing rate constants based
    % on the total edit statistic above
    k_total_edit_old = m1.Parameters(8).Value + m1.Parameters(9).Value;
    k_modifier = k_total_edit/k_total_edit_old;
    m1.Parameters(8).Value = m1.Parameters(8).Value * k_modifier;
    m1.Parameters(9).Value = m1.Parameters(9).Value * k_modifier;
    
    % Run the Model
    [Time, SimData, Names] = sbiosimulate(m1, csObj, doseObj);
    
    % Capture RNP
    RNP = SimData(:,7);
    
    % Calculate the Flux by multiplying rates of the 4 editing equations
    k_total = m1.Parameters(8).Value + m1.Parameters(9).Value;
    Cells_Total = SimData(:,1) + SimData(:,2);
    dE_dt = k_total .* Cells_Total .* SimData(:,7);
    Total_Cells = SimData(:,1) + SimData(:,2) + SimData(:,3) + SimData(:,4) + SimData(:,5) + SimData(:,6);
    dE_dt_Norm = dE_dt ./ Total_Cells;
end