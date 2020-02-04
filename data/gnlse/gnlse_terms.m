%% Compute terms in GNLSE (nondimensionalized input)
% Jared Callaham, 2019

clear all; close all; clc;
load('gnlse_nondim.mat')

%%
order=10; % Highest derivative to include
dt = t(2) - t(1);
dx = x(2) - x(1);

it = @(field) 10*log10(abs(field).^2); % log scale intensity

n = length(t); % grid size
freq = fftshift(2*pi*(-n/2:n/2-1)'/(n*dt))'; % frequency grid (centered)

uhat = fftshift( ifft( u, n, 2), 2 )/dt;


%% Plot energy content in spectral domain
V = 2*pi*(-n/2:n/2-1)'/(n*dt); % frequency grid

%% Travers code for RHS computation (modified to nondimensionalized)
B = 0;
for i = 1:length(alphas)        % Taylor expansion of betas
  B = B + alphas(i)/factorial(i+1).*V.^(i+1);
end
L = 1i*B;            % linear operator
W = V + w0;                    % for shock W is true freq (center-frequency shifted)
rw = n*ifft(fftshift(rt.'));   % frequency domain Raman term
L = fftshift(L).';
W = fftshift(W).'; % shift to fft space

rhs_full = zeros(size(uhat));
for k=1:length(x)
    % Invert variable transform
    fieldW = fftshift(uhat(k,:))*dt;  % scale
    fieldT = fft(fieldW);         % time domain field
    
    IT = abs(fieldT).^2;                % time domain intensity
    RS = dt*fr*fft(ifft(IT).*rw.');       % Raman convolution (temporal domain?)
    M = ifft(fieldT.*((1-fr).*IT + RS));  % full response function
    rhs_full(k, :) = 1i*W.*M/w0;   % full RHS of Eq. (3.13)

    % Transform back
    rhs_full(k,:) = fftshift( L.*fieldW + rhs_full(k, :) )./dt;
end


%% All RHS terms
ux = fft( dt*fftshift(rhs_full, 2), n, 2);


%% Split up linear operator (FFTSHIFTED VS ABOVE)
L_split = zeros(order, n);
for k = 2:order       % Taylor expansion of betas
  L_split(k-1, :) = 1i*alphas(k-1)/factorial(k).*V.^(k);
end

duhat_dt = {};
du_dt = {};
rhs_split = zeros(size(uhat));
for k=2:order
    duhat_dt{k} = zeros(size(uhat));
    for j=1:length(x)
        duhat_dt{k}(j, :) = L_split(k-1, :).*uhat(j, :);
    end
    rhs_split = rhs_split + duhat_dt{k};
    du_dt{k} = fft( dt*fftshift(duhat_dt{k}, 2), n, 2);
end


%% Nonlinear part of RHS
rhs_nonlinear = zeros(size(uhat));
for k=1:length(x)
    
    % Invert variable transform
    fieldW = fftshift(uhat(k,:))*dt;  % scale
    fieldT = fft(fieldW);         % time domain field
    
    IT = abs(fieldT).^2;                % time domain intensity
    RS = dt*fr*fft(ifft(IT).*rw.');       % Raman convolution (temporal domain?)
    M = ifft(fieldT.*((1-fr).*IT + RS));  % full response function
    rhs_full(k, :) = 1i*W.*M/w0;   % full RHS of Eq. (3.13)

    % Transform back
    rhs_nonlinear(k,:) = fftshift( rhs_nonlinear(k,:) )./dt;
end

kerr = 1i*u.*abs(u).^2;   % Cubic (kerr) nonlinearity
raman = fft( dt*fftshift(rhs_nonlinear, 2), n, 2) - kerr;  % Raman term (without Kerr)


%% NEED TO CLEAN UP... just save for now
save('../gnlse_nondim.mat', 'u', 'ux', 'du_dt', 'kerr', 'raman', 'x', 't', 'alphas')