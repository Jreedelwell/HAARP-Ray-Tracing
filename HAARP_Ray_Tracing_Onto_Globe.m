%
% Purpose :
%   Example of using raytrace_3d for a fan of rays.
%
% Calling sequence :
%   ray_test_3d
%
% Inputs :
%   None
%
% Outputs :
%   None
%
% Modification History:
%   V1.0  M.A. Cervera  07/12/2009
%     Initial version.
%
%   V1.1  M.A. Cervera  12/05/2009
%     Uses 'parfor' to parallelize the computation if the parallel computing
%     tool box is available
%
%   V1.3  M.A. Cervera  19/05/2011
%     More efficient handling of ionospheric  and geomagnetic grids grids in
%     call to raytrace_3d
%
%   V2.0 M.A. Cervera  03/05/2016
%     Modified to make use of multi-threaded raytrace_3d. IRI2016 is now used
%     to generate the ionosphere.
%

%
% setup general stuff
%
UT = [2020 10 15 18 30];           % UT - year, month, day, hour, minute
speed_of_light = 2.99792458e8;
R12 = 100;
elevs = [1:1:180]; % initial elevation of rays
freq = 9.6;
freqs = freq.*ones(size(elevs));   % frequency (MHz)
ray_bears = azimuth(62.39232,-145.14217,70.47388,-68.58611)*ones(size(elevs)); % initial bearing of rays

%New Changes
brngs = unique(ray_bears)-30:5:unique(ray_bears) + 30;
[brng_2d, elv_2d] = meshgrid(elevs, brngs);
brng_list = transpose(brng_2d(:));
elv_list = transpose(elv_2d(:));
freq_list=freq*ones(size(brng_list));
%

origin_lat = 62.39232;             % latitude of the start point of rays
origin_long = -145.14217;            % longitude of the start point of rays
origin_ht = 0.0;                % altitude of the start point of rays
doppler_flag = 1;               % interested in Doppler shift

fprintf( ['\n' ...
    'Example of 3D magneto-ionic numerical raytracing for a WGS84 ellipsoidal' ...
    ' Earth\n\n'])

%
% generate ionospheric, geomagnetic and irregularity grids
%
ht_start = 60;          % start height for ionospheric grid (km)
ht_inc = 2;             % height increment (km)
num_ht = 201;
lat_start = 57.39232;
lat_inc = 0.3;
num_lat = 101.0;
lon_start= -150.15;
lon_inc = 1.0;
num_lon = 101.0; %CANNOT GO HIGHER FOR SOME REASON
iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht, ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc = 10;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1.0;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);
geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];


tic
fprintf('Generating ionospheric and geomag grids... ')
[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
    gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
    geomag_grid_parms, doppler_flag);
toc

%FOF2 modification
iono_pf_grid  = iono_pf_grid / 2.1;
iono_pf_grid_5  = iono_pf_grid_5 / 2.1;

fprintf('\n')

% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;

%
% call raytrace
%
nhops = 4;                  % number of hops
tol = [1e-7 0.01 25];       % ODE solver tolerance and min max stepsizes
num_elevs = length(elv_list);


% Generate the O mode rays
OX_mode = 1;

fprintf('Generating %d O-mode rays ...', num_elevs);
tic
[ray_data_O, ray_O, ray_state_vec_O] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elv_list, brng_list, freq_list, ...
    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
    collision_freq, iono_grid_parms, Bx, By, Bz, ...
    geomag_grid_parms);

NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
    [ray_data_O.NRT_elapsed_time], NRT_total_time)


for rayId=1:num_elevs
    num = length(ray_O(rayId).lat);
    ground_range = zeros(1, num);
    lat = ray_O(rayId).lat;
    lon = ray_O(rayId).lon;
    ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
        origin_long,'wgs84')/1000.0;
    ray_O(rayId).ground_range = ground_range;
end


% Generate the X mode rays - note in the raytrace_3d call the ionosphere does
% not need to be passed in again as it is already in memory
OX_mode = -1;

fprintf('Generating %d X-mode rays ...', num_elevs);
tic
[ray_data_X, ray_X, ray_sv_X] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elv_list, brng_list, freq_list, ...
    OX_mode, nhops, tol);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
    [ray_data_X.NRT_elapsed_time], NRT_total_time)

for rayId=1:num_elevs
    num = length(ray_X(rayId).lat);
    ground_range = zeros(1, num);
    lat = ray_X(rayId).lat;
    lon = ray_X(rayId).lon;
    ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
        origin_long,'wgs84')/1000.0;
    ray_X(rayId).ground_range = ground_range;
end


% Generate the rays for the case where the magnetic field is ignored  - note
% in the raytrace_3d call the ionosphere does not need to be passed in again
% as it is already in memory
OX_mode = 0;

fprintf('Generating %d ''no-field'' rays ...', num_elevs);
tic
[ray_data_N, ray_N, ray_sv_N] = ...
    raytrace_3d(origin_lat, origin_long, origin_ht, elv_list, brng_list, freq_list, ...
    OX_mode, nhops, tol);
NRT_total_time = toc;
fprintf('\n   NRT-only execution time = %f, Total mex execution time = %f\n\n', ...
    [ray_data_N.NRT_elapsed_time], NRT_total_time)

for rayId=1:num_elevs
    num = length(ray_N(rayId).lat);
    ground_range = zeros(1, num);
    lat = ray_N(rayId).lat;
    lon = ray_N(rayId).lon;
    ground_range(2:num) = latlon2raz(lat(2:num), lon(2:num), origin_lat, ...
        origin_long,'wgs84')/1000.0;
    ray_N(rayId).ground_range = ground_range;
end

fprintf('\n')

% finished ray tracing with this ionosphere so clear it out of memory
clear raytrace_3d

%plotting the rays using 3D Earth Example
%Textured 3D Earth example
%
% Ryan Gray
% 8 Sep 2004
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013

%Options

space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
GMST0 = []; % Don't set up rotatable globe (ECEF)
%GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.

image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

% Mean spherical earth

erad    = 6371008.7714; % equatorial radius (meters)
prad    = 6371008.7714; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%Create figure

figure('Color', space_color);

hold on;

% Turn off the normal axes

set(gca, 'NextPlot','add', 'Visible','off');

axis equal;
axis auto;

% Set initial view

view(0,30);

axis vis3d;

%Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);


globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);


if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%Texturemap the globe

% Load Earth image for texture map

cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');


%Plotting the rays in O-mode
for ii = 1:num_elevs

    ray_Lon = ray_O(ii).lon;
    ray_Lat = ray_O(ii).lat;
    ray_Ht = ray_O(ii).height+200;

    [f,g,h] = geodetic2ecef(wgs84Ellipsoid,ray_Lat,ray_Lon,ray_Ht);
    plot3(f,g,h, 'b.', ...
        'markersize', 5);
 

end

%Ploting the points perpindicular to B-field
for rayId=1:num_elevs
    perp_B_idx = ray_O(rayId).wavenorm_B_angle > 89 & ray_O(rayId).wavenorm_B_angle < 91;
    if sum(perp_B_idx) ~= 0
        disp('found some points')
        lats = ray_O(rayId).lat(perp_B_idx);
        lons = ray_O(rayId).lon(perp_B_idx);
        hts = ray_O(rayId).height(perp_B_idx)+200;

        [w,e,r] = geodetic2ecef(wgs84Ellipsoid,lats,lons,hts);

        plot3(w,e,r, 'r.', 'markersize', 15);

    end
end


hold on