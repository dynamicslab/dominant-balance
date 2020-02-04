%% Nondimensionalize output of GNLSE simulation
%
% Jared Callaham, 2018

clear all; close all; clc;
load('gnlse.mat')

% Convert everything to nm and ps
betas = betas*1e-9;  % ps/m to ps/nm
Z = 1e9*Z;            % m to nm
gamma = 1e-9*gamma;   % 1/W/m to 1/W/nm

%% Define scaling constants

% Natural units
% T0 = 1/w0;
% Z0 = wavelength;
% A0 = sqrt(gamma*wavelength);

% Unity scaling
% T0 = t0;
% Z0 = 1e9*flength;  % Including m to nm
% A0 = sqrt(power);

% Soliton scaling
A0 = sqrt(power);
Z0 = 1/((1-fr)*gamma*A0^2);
T0 = sqrt(Z0*abs(betas(1)));

%% Space/time variables
t = T/T0;
x = Z/Z0;
u = AT/A0;
w0 = w0*T0;  % Nondimensional frequency

%% Taylor coefficients
alphas = zeros(size(betas));

% Only nondimensionalization
for k=1:length(alphas)
	alphas(k) = (Z0*betas(k))/(T0^(k+1));
end

%% Raman kernel
rt = gamma*T0*A0^2*Z0*(tau1^2+tau2^2)/tau1/tau2^2*exp(-T/tau2).*sin(T/tau1);
rt(t<0) = 0;          % Heaviside step function

%% Plot field
it = @(field) 10*log10(abs(field).^2); % log scale intensity

figure()
lIT = it(u); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(t, x, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
xlim([-100, 1400]);
shading interp;
xlabel('Delay'); ylabel('Distance');
colorbar()

%% Save output

save('gnlse_nondim.mat', 'x', 't', 'u', 'c', 'w0', 'rt', 'alphas',...
    'tau1', 'tau2', 'fr', 'gamma', 'wavelength', 'A0', 'T0', 'Z0')