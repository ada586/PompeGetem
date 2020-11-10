clear all
%% Script to evaluate the normalized dEdt for all the k values generated in the k_Generator script

% Load Structure with the models and the k values
load('ValidationdEdtData\ValidationStructureFile.mat');

% Calulate what N_bootstrap is
N_Bootstrap = length(ValidationModel(1).BootstrappedReplicates);

for i=1:4
    sbioloadproject(ValidationModel(i).FileNames);
    sbioaccelerate(m1);
    tic;
    for j=1:N_Bootstrap
        [Time, RNP, dE_dt_Norm] = RNP_Flux_Dynamics(m1, ValidationModel(i).BootStrappedk(j));
        ValidationModel(i).NormalizeddEdt(j) = trapz(Time, dE_dt_Norm)/max(Time);
        ValidationModel(i).SimulationData(j).Time = Time;
        ValidationModel(i).SimulationData(j).RNP = RNP;
        ValidationModel(i).SimulationData(j).dE_dt_Norm = dE_dt_Norm;
        [i,j]
    end
end

save('ValidationdEdtData/ValidationStructureFile_dEdt_Norm.mat','ValidationModel')