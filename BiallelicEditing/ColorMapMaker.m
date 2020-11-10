function Colormaps = ColorMapMaker(Granularity)
ERT_Healing = [0.3267 0.0958 0.5387];
N_CMapColors = Granularity * 100;
C_Boundary = round(ERT_Healing(1) * N_CMapColors);
Colormaps.C.Array = zeros(N_CMapColors, 3);
for i=1:C_Boundary
    Colormaps.C.Array(i,:) = [1 i/C_Boundary i/C_Boundary];
end
for i = C_Boundary:N_CMapColors
    Colormaps.C.Array(i,:) = [(N_CMapColors - i)/(N_CMapColors - C_Boundary) (N_CMapColors - i)/(N_CMapColors - C_Boundary) 1];
end

S_Boundary = round(ERT_Healing(2) * N_CMapColors);
Colormaps.S.Array = zeros(N_CMapColors, 3);
for i=1:S_Boundary
    Colormaps.S.Array(i,:) = [1 i/S_Boundary i/S_Boundary];
end
for i = S_Boundary:N_CMapColors
    Colormaps.S.Array(i,:) = [(N_CMapColors - i)/(N_CMapColors - S_Boundary) (N_CMapColors - i)/(N_CMapColors - S_Boundary) 1];
end

L_Boundary = round(ERT_Healing(3) * N_CMapColors);
Colormaps.L.Array = zeros(N_CMapColors, 3);
for i=1:L_Boundary
    Colormaps.L.Array(i,:) = [1 i/L_Boundary i/L_Boundary];
end
for i = L_Boundary:N_CMapColors
    Colormaps.L.Array(i,:) = [(N_CMapColors - i)/(N_CMapColors - L_Boundary) (N_CMapColors - i)/(N_CMapColors - L_Boundary) 1];
end
end