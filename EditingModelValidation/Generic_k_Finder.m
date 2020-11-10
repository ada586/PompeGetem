%% Function that returns a k value given a model and the percentage of edited cells as the input
function [k_total_edit] = Generic_k_Finder(m1, PercentEdit)
    % Start by initiating the model as written
    csObj = getconfigset(m1, 'Active');
    csObj.TimeUnits = 'day';
    doseObj = m1.getdose;
    
    k_total_edit = m1.Parameters(8).Value + m1.Parameters(9).Value;
    % Start 100 fold iterative loop of taking in k, finding PercentEdit and then modifying k till answer
    for i = 1:100
    k_total_edit_old = m1.Parameters(8).Value + m1.Parameters(9).Value;
    k_multiplier = k_total_edit/k_total_edit_old;
    m1.Parameters(8).Value = m1.Parameters(8).Value * k_multiplier;
    m1.Parameters(9).Value = m1.Parameters(9).Value * k_multiplier;
    
    [Time SimData Names] = sbiosimulate(m1, csObj, doseObj);
    Total_Cells = SimData(:,1) + SimData(:,2) + SimData(:,3) + SimData(:,4) + SimData(:,5) + SimData(:,6);
    Precise_Edited_Cells = SimData(:,3) + SimData(:,4);
    n_datapoints = length(Time);
    PercentEditEstimate = Precise_Edited_Cells(n_datapoints) / Total_Cells(n_datapoints) * 100;
    k_Modifier(i) = PercentEdit/PercentEditEstimate;
    k_total_edit = k_total_edit * k_Modifier(i);
    end
end