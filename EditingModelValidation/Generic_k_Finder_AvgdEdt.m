%% Function that returns a k value given an accelerated model and the percentage of edited cells as the input
function [k_total_edit k_Modifier k_History] = Generic_k_Finder_AvgdEdt(m1, AvgdEdt)
    % Start by initiating the model as written
    csObj = getconfigset(m1, 'Active');
    csObj.TimeUnits = 'day';
    doseObj = m1.getdose;
    
    k_total_edit = m1.Parameters(8).Value + m1.Parameters(9).Value;
    % Start 100 fold iterative loop of taking in k, finding PercentEdit and then modifying k till answer
    for i = 1:25
        k_total_edit_old = m1.Parameters(8).Value + m1.Parameters(9).Value;
        k_multiplier = k_total_edit/k_total_edit_old;
        m1.Parameters(8).Value = m1.Parameters(8).Value * k_multiplier;
        m1.Parameters(9).Value = m1.Parameters(9).Value * k_multiplier;
        
        [Time, RNP, dE_dt_Norm] = RNP_Flux_Dynamics(m1, k_total_edit);
        dE_dt_Norm_Avg = trapz(Time, dE_dt_Norm)/max(Time);
        k_Modifier(i) = AvgdEdt/dE_dt_Norm_Avg
        k_total_edit = k_total_edit * k_Modifier(i);
        k_History(i) = k_total_edit;
    end
end