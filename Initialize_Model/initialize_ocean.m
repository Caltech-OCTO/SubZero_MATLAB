function [ocean, heat_flux, h0]=initialize_ocean(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.U = 0.5;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=2e3; % in meters

Lx = 1e5 / 2; kx = pi/Lx; Ly = 1e5 / 2; ky = pi/Ly;
Xo=-Lx:dXo:Lx; Yo=-Ly:dXo:Ly; 

%calculating ocean velocity field
% Assume Yo goes from -Ly to Ly, Nx = length(Xo), Ny = length(Yo)
%% Original domain
Ny = length(Yo);
Nx = length(Xo);
half_ny = floor(Ny/2);

% Original piecewise linear profile along y
u_vec = [ ...
    linspace(0, ocean.U, half_ny)'; ...            % from y_min to y=0
    ocean.U; ...                                   % center point y=0
    linspace(ocean.U, 0, Ny - half_ny - 1)' ...   % from y=0 to y_max
];

%% Original spacing
dy = Yo(2) - Yo(1);
dx = Xo(2) - Xo(1);

%% Extended coordinates for ±1.5Ly and ±1.5Lx
Yo_ext = (-2*Ly : dy : 2*Ly).';
Xo_ext = (-2*Lx : dx : 2*Lx);
[Xocn_ext, Yocn_ext]=meshgrid(Xo_ext,Yo_ext);

Ny_ext = length(Yo_ext);
Nx_ext = length(Xo_ext);

%% Map extended y to periodic indices into u_vec
y_offset = Yo_ext - 0;  % relative to center at y=0

% Floating point index in u_vec
y_index_float = (y_offset + Ly)/dy + 1;

% Wrap indices modulo Ny to repeat pattern
y_index_mod = mod(y_index_float - 1, Ny) + 1;  % range [1, Ny]

% Round to nearest integer and clip to valid indices
y_indices = round(y_index_mod);
y_indices(y_indices < 1) = 1;
y_indices(y_indices > Ny) = Ny;

% Build extended profile
u_vec_ext = u_vec(y_indices);



%% Repeat along x
Uocn = repmat(u_vec_ext, 1, Nx_ext);
Vocn = zeros(size(Uocn));

ocean.Xo=Xo_ext;
ocean.Yo=Yo_ext;
ocean.Xocn = Xocn_ext;
ocean.Yocn = Yocn_ext;
ocean.kx = kx;
ocean.ky = ky;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;

%Calculate heat flux and how much ice would grow between creation of new
%floes
k = 2.14; %Watts/(meters Kelvin)
dt = 20; %seconds
Ta = -1; %Celcius
To = 0; %Celcius
rho_ice = 920; %kg/m^3
L = 2.93e5; % Joules/kg

h0 = real(sqrt(2*k*dt*nDTOut*(To-Ta)/(rho_ice*L)));
heat_flux = k*(Ta-To)/(rho_ice*L); 
h0 = mean(h0(:)); %Thickness of newly created sea ice;

end
