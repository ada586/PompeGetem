clear all
%% Script to make figures for the dEdt for the training set with bootstrapped replicates
load('ValidationdEdtData\ValidationStructureFile_dEdt_Norm.mat')

for i=1:4
    sbioloadproject(ValidationModel(i).FileNames);
    csObj = getconfigset(m1);
    csObj(1).StopTime = 35;
    DoseObj = getdose(m1);
    k_total_edit(i) = mean(ValidationModel(i).BootStrappedk);
    [Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(i));
    figure;
    plot(Time, dEdt_Norm);
    set(gca,'YLim', [0 .120]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end