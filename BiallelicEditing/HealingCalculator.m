function [Time Cardiac Skeletal Liver] = HealingCalculator(m1, csObj, DoseObj)
    % Function to run the model and perform the calculations as needed, to be
    % called in loops in heatmaps
    [Time SimData Names] = sbiosimulate(m1, csObj, DoseObj);
    % Calculate Cardiac healing
    Cardiac = SimData(:,34)./(SimData(:,33) + SimData(:,34));
    Skeletal = SimData(:,36)./(SimData(:,35) + SimData(:,36));
    Liver = SimData(:,37);
end