clear all;

Granularity = 101;
Progenitor_Dose = 100e6 * linspace(2, 0, Granularity);
Half_Life = linspace(50, 300, Granularity);
GAA_Boost = linspace(0, 1.5, Granularity);

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Functionalized(Progenitor_Dose(i), 0.75, Half_Life(j)); 
    end
    i
end

Endpoint_Dose_v_Half_Life = Endpoint;

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Functionalized(Progenitor_Dose(i), GAA_Boost(j), 250); 
    end
    i
end

Endpoint_Dose_v_GAA_Boost = Endpoint;

for i=1:Granularity
    for j=1:Granularity
    Endpoint(i,j,:) = Cell_Therapy_Pompe_Disease_Graph_Generator_Functionalized(100e6, GAA_Boost(Granularity + 1 - i), Half_Life(j)); 
    end
    i
end

Endpoint_GAA_Boost_v_Half_Life = Endpoint;

save('Cell_Therapy_Heatmap_data_101','Endpoint_Dose_v_GAA_Boost', 'Endpoint_Dose_v_Half_Life', 'Endpoint_GAA_Boost_v_Half_Life', 'Progenitor_Dose', 'GAA_Boost', 'Half_Life');