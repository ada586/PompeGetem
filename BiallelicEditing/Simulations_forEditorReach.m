clear all
%% Script to run experiments required for reviewer 4
% Reach of the genome editors
%% Load project, configset and doses
sbioloadproject('SimBiologyProjectFiles\Pompe_Model_v28.sbproj')
csObj = getconfigset(m1);
DoseObj = getdose(m1);
%% Run the model for 29 days, with one editor only
csObj(1).StopTime = 29;
[T D N] = sbiosimulate(m1, csObj(1), DoseObj(2));

Unedited_Indices = [1 2 19 20];
Allele_1_Indices = [3 4 7 8 21 22];

TotalCellArray = zeros(length(T),1);
for j=1:26
    TotalCellArray = TotalCellArray + D(:,j);
end
figure;

% Rate of liver genome editing
k_total_edit = m1.Parameters(15).Value;
Edit_Flux = k_total_edit * D(:,27) .* (D(:,1)+D(:,2)+D(:,19)+D(:,20));
Norm_Edit_Flux = Edit_Flux ./ TotalCellArray;
p = plot(T, Norm_Edit_Flux * 100); hold on;
CellsReached = trapz(T, Norm_Edit_Flux * 100)/T(end);
% Precision = m1.Parameters(13).Value/(1+m1.Parameters(13).Value);
% k_precise = Precision * k_total_edit;
% Prec_Edit_Flux = k_precise * D(:,27) .* (D(:,1)+D(:,2)+D(:,19)+D(:,20));
% Prec_Norm_Edit_Flux = Prec_Edit_Flux ./ TotalCellArray;
% plot(T, Prec_Norm_Edit_Flux * 100); hold all;

% Rate of liver genome change
n = length(T);
Edited = zeros(length(T),1);
for j=Allele_1_Indices
    Edited = Edited + D(:,j);
end

for i=2:n-2
    DerivEdited(i-1) = (Edited(i+2) - Edited(i))/(T(i+2) - T(i));
end
TotalCellArray = TotalCellArray';
Norm_DerivEdited = DerivEdited ./ TotalCellArray(2:n-2);
Treshape = T(2:n-2);
p = plot(Treshape, Norm_DerivEdited * 100); hold on;
CellsEdited = trapz(Treshape,Norm_DerivEdited * 100)/T(n-1);

CC_Indices = [19 20 21 22];
CCed = zeros(length(T),1);
for j=CC_Indices
    CCed = CCed + D(:,j);
end

for i=2:n-2
    DerivCCed(i-1) = (CCed(i+2) - CCed(i))/(T(i+2) - T(i));
end
Norm_DerivCCed = DerivCCed ./ TotalCellArray(2:n-2);
Treshape = T(2:n-2);
p = plot(Treshape, Norm_DerivCCed * 100);
CellsEdited = trapz(Treshape,Norm_DerivCCed * 100)/T(n-1);

set(gca,'YLim',[0 1.5]);
set(gca,'XLim',[0 30]);
figure;

%% Percentage Edited vs cross corrected
Precisely_Edited = [3 4];
Imprecisely_Edited = [7 8 21 22];
CrossCorrected = [19 20 21 22];

Precise = zeros(length(T),1);
for j = Precisely_Edited
    Precise = Precise + D(:,j); 
end

Imprecise = zeros(length(T),1);
for j = Imprecisely_Edited
    Imprecise = Imprecise + D(:,j); 
end

CC = zeros(length(T),1);
for j = CrossCorrected
    CC = CC + D(:,j); 
end

Norm_Precise = Precise ./ TotalCellArray';
Norm_Imprecise = Imprecise ./ TotalCellArray';
Norm_CC = CC ./ TotalCellArray';

plot(T, Norm_Precise, T, Norm_Imprecise, T, Norm_CC);
