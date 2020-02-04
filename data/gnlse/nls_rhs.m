%% Check computation of RHS terms vs Dudley code
% Jared Callaham, 2018

clear all; close all; clc;
addpath(genpath('../utils'));
load('../data/gnlse.mat')

%%


order=10; % Highest derivative to include
dT = T(2) - T(1);
dZ = Z(2) - Z(1);

it = @(field) 10*log10(abs(field).^2); % log scale intensity

n = length(T); % grid size
freq = fftshift(2*pi*(-n/2:n/2-1)'/(n*dT))'; % frequency grid (centered and nondimensionalized)
 
%% Plot energy content in temporal domain (from Travers code)
figure();

lIT = it(AT); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(T, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Delay'); ylabel('Distance');
%xlim(xrange);
colorbar()

%% Plot energy content in spectral domain
V = 2*pi*(-n/2:n/2-1)'/(n*dT); % frequency grid

figure()
lIT = it(AW); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
%xlim(xrange);
colorbar()

%% Travers code for RHS computation
B = 0;
for i = 1:length(betas)        % Taylor expansion of betas
  B = B + betas(i)/factorial(i+1).*V.^(i+1);
end
L = 1i*B;            % linear operator
W = V + w0;                % for shock W is true freq
RW = n*ifft(fftshift(RT.'));   % frequency domain Raman term
L = fftshift(L).';
W = fftshift(W).'; % shift to fft space

%% === 
rhs_full = zeros(size(AW));
for k=1:length(Z)
    % Invert variable transform
    fieldW = fftshift(AW(k,:))*dT;  % scale
    fieldT = fft(fieldW);         % time domain field
    
    IT = abs(fieldT).^2;                % time domain intensity
    RS = dT*fr*fft(ifft(IT).*RW.');       % Raman convolution (temporal domain?)
    M = ifft(fieldT.*((1-fr).*IT + RS));  % response function
    rhs_full(k, :) = 1i*(gamma/w0)*W.*M;   % full RHS of Eq. (3.13)

    % Transform back
    rhs_full(k,:) = fftshift( L.*fieldW + rhs_full(k, :) )./dT;
end

rhs_linear = zeros(size(AW));
for k=1:length(Z)
    rhs_linear(k,:) = fftshift(L).*AW(k,:);
end

rhs_nonlinear = zeros(size(AW));
for k=1:length(Z)
    % Invert variable transform
    fieldW = fftshift(AW(k,:))*dT;  % scale
    fieldT = fft(fieldW);         % time domain field
    
    IT = abs(fieldT).^2;                % time domain intensity
    RS = dT*fr*fft(ifft(IT).*RW.');       % Raman convolution (temporal domain?)
    M = ifft(fieldT.*((1-fr).*IT + RS));  % response function
    rhs_nonlinear(k, :) = 1i*(gamma/w0)*W.*M;   % full RHS of Eq. (3.13)

    % Transform back
    rhs_nonlinear(k,:) = fftshift( rhs_nonlinear(k,:) )./dT;
end

%% Finite-difference spatial derivative estimate (LHS of equation) - spectral domain
rhs_fd = zeros(size(AW));
% Forward difference at start
rhs_fd(1, :) = ( AW(2, :) - AW(1, :) )/dZ;
% Backwards difference at end
rhs_fd(end, :) = ( AW(end, :) - AW(end-1, :) )/dZ;

% Central difference in between
for j=2:length(Z)-1
    rhs_fd(j, :) = ( AW(j+1, :) - AW(j-1, :) )/(2*dZ);
end

figure()
subplot(121)
lIT = it(rhs_full); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
%xlim(xrange);
title('Dudley code')
colorbar()

subplot(122)
lIT = it(rhs_fd); % log scale temporal intensity
%mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
%xlim(xrange);
title('Finite-difference')
colorbar()

%% Try to break off derivative terms

% Split up linear operator (FFTSHIFTED VS ABOVE)
L_split = zeros(order, n);
for k = 2:order        % Taylor expansion of betas
  L_split(k-1, :) = 1i*betas(k-1)/factorial(k).*V.^(k);
end

AWz = {};
rhs_split = zeros(size(AW));
for k=2:order
    AWz{k} = zeros(size(AW));
    for j=1:length(Z)
        AWz{k}(j, :) = L_split(k-1, :).*AW(j, :);
    end
    rhs_split = rhs_split + AWz{k};
end

% Compare all terms in spectral domain
figure()
subplot(221)
lIT = it(rhs_full); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
title('Full Dudley code')
colorbar()

subplot(222)
lIT = it(rhs_nonlinear); % log scale temporal intensity
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
title('Dudley nonlinear terms')
colorbar()

subplot(223)
lIT = it(rhs_split); % log scale temporal intensity
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
title('Sum of separate linear terms')
colorbar()

subplot(224)
lIT = it(rhs_linear); % log scale temporal intensity
pcolor(V, Z, lIT);          % plot as pseudocolorap
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Frequency'); ylabel('Distance');
title('Dudley linear terms')
colorbar()


%% Compare derivatives in time domain
figure()
subplot(221)
field = fft( dT*fftshift(rhs_full, 2), n, 2);
lIT = it(field); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
pcolor(T, Z, lIT);          % plot as pseudocolor map
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Delay'); ylabel('Distance');
title('Full Dudley code')
colorbar()

subplot(222)
field = fft( dT*fftshift(rhs_nonlinear, 2), n, 2);
lIT = it(field); % log scale temporal intensity
pcolor(T, Z, lIT);          % plot as pseudocolor map
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Delay'); ylabel('Distance');
title('Dudley nonlinear terms')
colorbar()

subplot(223)
field = fft( dT*fftshift(rhs_split, 2), n, 2);
lIT = it(field); % log scale temporal intensity
pcolor(T, Z, lIT);          % plot as pseudocolor map
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Delay'); ylabel('Distance');
title('Sum of separate linear terms')
colorbar()

subplot(224)
field = fft( dT*fftshift(rhs_linear, 2), n, 2);
lIT = it(field); % log scale temporal intensity
pcolor(T, Z, lIT);          % plot as pseudocolor map
caxis([mlIT-40.0, mlIT]);
shading interp;
xlabel('Delay'); ylabel('Distance');
title('Dudley linear terms')
colorbar()

