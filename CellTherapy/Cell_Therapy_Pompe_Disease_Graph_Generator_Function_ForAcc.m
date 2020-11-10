function Endpoint = Cell_Therapy_Pompe_Disease_Graph_Generator_Function_ForAcc(CCModel, Progenitor_Dose, GAA_Boost, Half_Life)
%% Retrieve Current dosage information
All_Dose = getdose(CCModel);
CCModel_dose = [All_Dose(5)] ;
csObj = getconfigset(CCModel);

%% Change parameters
All_Dose(5).Amount = Progenitor_Dose;
CCModel.Parameters(18).Value = GAA_Boost;
CCModel.Parameters(27).Value = 2 * log(2)/Half_Life - CCModel.Parameters(26).Value;

%% Run Model
CCModeldata = sbiosimulate(CCModel, csObj(1), CCModel_dose);

%% Plotting the percentage of healthy cells
Time = CCModeldata.Time;
Heart_Total = CCModeldata.Data(:,33) + CCModeldata.Data(:,34); 
HealthyCardiac = CCModeldata.Data(:,34)./Heart_Total*100;

Skeletal_Total = CCModeldata.Data(:,35) + CCModeldata.Data(:,36); 
HealthySkeletal = CCModeldata.Data(:,36)./Skeletal_Total*100;

HealthyLiver = CCModeldata.Data(:,37) * 100;

Export_Data = [Time/30 HealthyCardiac HealthyLiver HealthySkeletal];
Export_Data_Usable = [(Time(50:end)/30 - 3) HealthyCardiac(50:end) HealthyLiver(50:end) HealthySkeletal(50:end)];

Endpoint = [HealthyCardiac(end) HealthySkeletal(end) HealthyLiver(end)];
end