clear all
%% Script to accelerate models and generate bootstrapped means and rate constants and write them to new file for all 4 validation models

% Adding in all the filenames into a structure
ValidationModel(1).FileNames = './SimbiologyProjectFiles/FAH_Editing_Project_v14_Yin_2014.sbproj';
ValidationModel(2).FileNames = './SimbiologyProjectFiles/FAH_Editing_Project_v16_Song_2018.sbproj';
ValidationModel(3).FileNames = './SimbiologyProjectFiles/FAH_Editing_Project_v16_Yin_2016.sbproj';
ValidationModel(4).FileNames = './SimbiologyProjectFiles/FAH_Editing_Project_v13_Yang_2016_AAV.sbproj';

% Adding in the editing data into the structure, setting seed of random
% number generator to 3 as it gives the closest range and mean for mRNA NP
% paper Yin 2016
rng(3);
ValidationModel(1).EditRange = [9.1 9.5];
ValidationModel(2).EditRange = [4.697674419, 5.953488372, 8.162790698, 11.62790698, 13.13953488, 13.93023256];
ValidationModel(3).EditRange = 1.1 * randn(4,1) + 6; % Figure 3g from Yin, 2016 paper, mean = 6, Stdev = 1.1, n=4
ValidationModel(4).EditRange = [5.87426141560978 7.041558505772704 8.366090632648522 9.880392593141984 19.15498529011842];

%% Generate Bootstrapped Replicates
N_Bootstrap = 10000;
for i=1:4
    ValidationModel(i).BootstrappedReplicates = bootstrp(N_Bootstrap, @mean, ValidationModel(i).EditRange);
end

%% Generate all the k values for each replicate
% Accelerate all of the Simbiology Projects and then calculate and store
% the rate constants
for i=1:4
    sbioloadproject(ValidationModel(i).FileNames);
    sbioaccelerate(m1);
    tic;
    for j=1:N_Bootstrap
        ValidationModel(i).BootStrappedk(j) = Generic_k_Finder(m1, ValidationModel(i).BootstrappedReplicates(j));
        [Time, RNP, dE_dt_Norm] = RNP_Flux_Dynamics(m1, ValidationModel(i).BootStrappedk(j));
        [i j]
    end
    ValidationModel(i).RunTime = toc;
end
save('ValidationdEdtData/ValidationStructureFile.mat','ValidationModel')