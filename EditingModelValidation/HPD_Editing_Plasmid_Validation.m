clear all
%% Validation of GETEM for HPD editing using nmCas9 Plasmid
% Load Validation Data structure and calculate the mouse RNP dose
load('ValidationdEdtData/ValidationStructureFile.mat')

% Calculate average and standard deviation of dEdt_Norm
k_mean = mean([ValidationModel(1).BootStrappedk, ValidationModel(2).BootStrappedk]);
k_min = min([ValidationModel(1).BootStrappedk, ValidationModel(2).BootStrappedk]);
k_max = max([ValidationModel(1).BootStrappedk, ValidationModel(2).BootStrappedk]);

% Load the PTEN editing with RNP
sbioloadproject('SimBiologyProjectFiles\HPD_Editing_Project_Ibrahim_2018.sbproj')

% k_total_edit = dEdt/RNP, calculate mean lower and upper bounds for
% k_total_edit

k_total_edit = [k_min k_mean k_max];

for i=1:3
    m1.Parameters(8).Value = k_total_edit(i);
    [Time, SimData, Names] = sbiosimulate(m1);
    N_Datapoints = length(Time);
    TotalCells(i) = SimData(N_Datapoints,1) + SimData(N_Datapoints,2) + SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditedCells(i) = SimData(N_Datapoints,3) + SimData(N_Datapoints,4) + SimData(N_Datapoints,5) + SimData(N_Datapoints,6);
    EditFraction(i) = EditedCells(i)/TotalCells(i);
end

save('ValidationdEdtData\HPD_Editing_Predictions.mat','EditFraction');

%% Graphing this 
csObj = getconfigset(m1);
DoseObj = getdose(m1);

csObj(1).StopTime = 49;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited)); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, Cells.Precise ./ (Cells.Precise + Cells.Imprecise + Cells.Unedited));
set(gca,'YLim', [0 .6]);
set(gca,'XLim', [0 49]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Ibraheim_Validation_TimeCourse.fig');
figure;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(1));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(3));
plot(Time, dEdt_Norm); hold all;
[Time, Cells, RNP, dEdt_Norm] = Indel_Dynamics(m1, csObj(1), DoseObj, k_total_edit(2));
plot(Time, dEdt_Norm);
set(gca,'YLim', [0 .120]);
set(gca,'XLim', [0 49]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
savefig('FigurePlots\Ibraheim_dEdt_Norm.fig');

PublishedData = [35 51 50 58 37];
