function dy = semireduced15neuron(t, y, Istim)
% In this version of the neuron, gNa = 0, all gating variables except those for INA and IK are set
% to their time->infty values, and Oc is set to its time->infty value

% time: msec
% current: nA
% potential: mV
% concentration: mM

%stimulating current
%Istim = 0 ;

% 
% V = -42;
% if t >= 5000 & t < 13000
%         V = 12;
% end

%V = -65;

%STATE VARIABLES:
V = y(1);  %membrane potential
Ca_i = y(2);  %intracellular Ca concentration, mM
Oc = y(3); %Ca buffer occupancy
m = y(4);
h = y(5);  %m, h: H-H gating variables for Na
% d = y(5);
% f = y(6);  %d, f: gating variables for fast Ca
% s = y(7);  %slow Ca
% b = y(9);  %non-specific cation
n = y(6);
l = y(7); %n, l: gating variables for delayed rectifier K

gamma = 0; % hack to fix Matlab bug where even though gamma is 
% defined as var in all_constants.m, it is still recognized as gamma.m
all_constants;


dy = zeros(size(y));

%~~~~~~~~~~Inward Currents~~~~~~~~~~~~~
%FAST SODIUM CURRENT
INa = gNa * m^3 * h * (V - ENa);
A_m = 0.08*(V+13.0) / (1 - exp((-V-13.0)/4.09));
A_h = 0.112 * exp((-30.0 - V)/10.0);
B_m = 2.15 * exp((-35.0 - V)/4.01);
B_h = 0.23 / (exp((9.65 - V)/23.9)+1);
tau_m = (A_m + B_m)^-1;
tau_h = (A_h + B_h)^-1;
dy(4) = (minf(V) - y(4))/tau_m;  %dm/dt
dy(5) = (hinf(V) - y(5))/tau_h;  %dh/dt



%FAST CALCIUM CURRENT
ICa = gCa * dinf(V)^2 * finf(V) * (V - ECa);
% A_d = 0.0063*(V + 10.81) / (1 - exp((-V-10.81)/5.03));
% A_f = 0.00325 * exp((10.0 - V)/7.57);
% B_d = 0.04 * exp((14.08 - V)/5.09);
% B_f = 0.029 / (1+ exp((20.29 - V)/5.4));
% tau_d = (A_d + B_d)^-1;
% tau_f = (A_f + B_f)^-1;
% dy(5) = (dinf(V) - y(5))/tau_d;  %dd/dt
% dy(6) = (finf(V) - y(6))/tau_f;  %df/dt

%SLOW INWARD CALCIUM CURRENT
% ISI = gSI * s * (V - ECa) * KMSI/(KMSI + Ca_i);
IISI = ISI(sinf(V), V, Ca_i);
% A_s = 0.0014 * (V - 54.0) / (1 - exp((-V+54.0)/12.63));
% B_s = 0.00013 * exp((-11.32-V)/16.8);
% tau_s = (A_s + B_s)^-1;
% dy(8) = (sinf(V) - y(8))/tau_s;  %ds/dt


%NON-SPECIFIC CATION CURRENT
INS = gNS * binf(V) * (V - ENS) * Ca_i/(Ca_i + KMNS);
% tau_b = 500 * (0.80 / (1 + exp((10.0 + V)/3.0)) + 0.20);
% dy(9) = (binf(V) - y(9))/tau_b;  %db/dt

%LEAKAGE CURRENT
IL = gL * (V - EL);


%~~~~~~~~~~Outward Currents~~~~~~~~~~~~~
%DELAYED RECTIFIER
IK = gK * n^4 * l * (V - EK);
A_n = 0.0035 * (V+17.0) / (1 - exp((-V-17.0)/3.0));
B_n = 0.04 * exp((-28.0 - V)/10.0);
tau_n = (A_n + B_n)^-1;
tau_l = 2000 * (0.90 / (1 + exp((28.0 + V)/3.0)) + 0.10);
tau_l = 6000 * (0.90 / (1 + exp((28.0 + V)/3.0)) + 0.10);
%tau_l = 15000;
dy(6) = (ninf(V) - y(6))/tau_n;  %dn/dt
dy(7) = (linf(V) - y(7))/tau_l;  %dl/dt

%ANOMALOUS RECTIFIER
Z = 1;  %Z was not defined in the paper; I assume it's a charge.
IR = gR * (V - EK + 5.66) / (1 + exp((V - EK - 15.3)*Z * F/(R*T)));


%~~~~~~~~~~Pumps and Exchangers~~~~~~~~~~~~~
ICaP = I_CaP * (Ca_i / (Ca_i + KMCaP));
DFin = (Na_i)^r * (Ca_o) * exp((r-2) * gamma * V *F/ (R*T));
DFout = (Na_o)^r * (Ca_i) * exp((r-2) * (gamma-1) * V * F/ (R*T));
S = 1 + DNaCa * (Ca_i*(Na_o)^r + Ca_o*(Na_i)^r);
INaCa = KNaCa * (DFin - DFout) / S;


%~~~~~~~~~~Internal Calcium Concentration~~~~~~~~~~~~~
% Oc = Oc_inf(Ca_i);
% dy(2) = (INaCa - IISI - ICa - ICaP - 0.197/gNS*(INS * (V - ECa)/(V - ENS))) / (2*Vol_i*F);  %dCa_i/dt
% dy(2) = dy(2) / (1 + nn * B_i * Oc^2 *kR / (kU*Ca_i^2));

%[INaCa/(2*Vol_i*F) -ISI/(2*Vol_i*F) -ICa/(2*Vol_i*F) -ICaP/(2*Vol_i*F) -nn*B_i*dy(3)]
% dy(2) = 380e-6 - y(2);

dy(3) = kU*y(2)*(1 - y(3)) - kR*y(3);  %dOc/dt
dy(2) = (INaCa - IISI - ICa - ICaP - 0.197/gNS*(INS * (V - ECa)/(V - ENS))) / (2*Vol_i*F) - nn * B_i * dy(3);  %dCa_i/dt

%~~~~~~~~~~Membrane Potential~~~~~~~~~~~~~
dy(1) = - (INa + ICa + IISI + INS + IK + IR + IL + INaCa + INaK + ICaP - Istim) / Cm;

%dy(1) = - (ICa + ISI + INS + IK + IR + IL + INaCa + INaK + ICaP - Istim) / Cm;
% dy(1) = -65 - y(1);
