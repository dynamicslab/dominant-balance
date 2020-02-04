% Download HYCOM Gulf of Mexico 2019 reanalysis from THREDDS database
% source: https://hycom.org/
% code: Jared Callaham (2018)
%
% Matlab NetCDF/OpenDAP docs:
%  https://www.mathworks.com/help/matlab/import_export/importing-network-common-data-form-netcdf-files-and-opendap-data.html
clear all; close all; clc;
 
%% Load expt_90.1 (1/1/2019-present hourly)
source = 'https://tds.hycom.org/thredds/dodsC/GOMu0.04/expt_90.1m000/data/hindcasts/2019?lat[0:1:345],lon[0:1:540],time[0:1:5063],water_u[0:1:5063][0][0:1:345][0:1:540],water_v[0:1:5063][0][0:1:345][0:1:540],water_temp[0:1:5063][0][0:1:345][0:1:540],salinity[0:1:5063][0][0:1:345][0:1:540],surf_el[0:1:5063][0:1:345][0:1:540]';
ncdisp(source)

%% Load expt_50.1 (1993-2012)
% Pull out Hurricane Rita path from September 2005
% URL: http://tds.hycom.org/thredds/dodsC/GOMu0.04/expt_50.1.html

source = 'http://tds.hycom.org/thredds/dodsC/GOMu0.04/expt_50.1?time[0:1:58410]';
t = ncread(source, 'time');
t_mask = (t > 49680) && (t < 50400); % September 2005

%%
lon = ncread(source, 'lon');  % 'x' coordinates on map
lat = ncread(source, 'lat');   % 'y' coordinates on map
t = ncread(source, 'time');
t = t(1:720); % First month of data

idx = 0;

%%
idx = 45;
for j=46:length(t)
    idx = idx + 1;
    fprintf('Loading %d/%d\n', j, length(t))
    uvel(:, :, idx) = (ncread(source, 'water_u', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % East-west velocity
    vvel(:, :, idx) = (ncread(source, 'water_v', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % North-south velocity
    temp(:, :, idx) = (ncread(source, 'water_temp', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % Water temperature
    sal(:, :, idx) = (ncread(source, 'salinity', [1, 1, 1, j], [Inf, Inf, 1, 1]))'; % Water salinity
    ssh(:, :, idx) = (ncread(source, 'surf_el', [1, 1, j], [Inf, Inf, 1]))'; % Sea surface height
end


%%
uvel = single(uvel);
vvel=single(vvel);
temp = single(temp);
sal = single(sal);
ssh = single(ssh);

[ny, nx, m] = size(uvel);
n = nx*ny;

%%
[X, Y] = meshgrid(lon, lat);

%% Convert to vorticity and save
vort = zeros(n, m);
for i=1:m
    [omega_z, ~]= curl(X,Y,uvel(:, :, i),vvel(:, :, i));
    vort(:, i) = reshape(omega_z(:, :, 1), [n, 1]);
end

save('./gom.mat', 'lat', 'lon', 'uvel', 'vvel', 'vort'...
    , 'temp', 'sal', 'ssh', 't')  % File will be ~1.1 Gb