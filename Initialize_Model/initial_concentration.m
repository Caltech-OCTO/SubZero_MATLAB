function [Floe, Nb] = initial_concentration(height, min_floe_size, Lx, Ly)
%% This function is used to generate the initial floe field

load("shear_flow_floes.mat");   % loads variable: polygons

numPolys = numel(polygons);
ps = polyshape.empty(numPolys, 0);

for i = 1:numPolys
    xy = polygons{i};
    xs = xy(:,1) - Lx;  % matlab -Lx to Lx, julia from 0 to 2Lx
    ys = xy(:,2) - Ly;  % matlab -Ly to Ly, julia from 0 to 2Ly
    ps(i) = polyshape(xs, ys);
end

Nb = 0;  % no boundaries in this simulation
Floe = [];

%Loop through all the polyshapes and create new floes
for ii = 1:numPolys
    floenew = initialize_floe_values(ps(ii), height);
    Floe = [Floe floenew];
end

areas = cat(1,Floe.area);
Floe(areas<min_floe_size)=[];
end

