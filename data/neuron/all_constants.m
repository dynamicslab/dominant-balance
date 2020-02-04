%CONSTANT PARAMETERS
gNa = 30; %uS
gCa = 10; %uS
gK  = 70; %uS
gNS = 0.5; %uS, non-specific

gSI = 0.747; %uS, slow inward Ca
gR  = 0.32; %uS, anomalous rectifier
gL  = 0.1; %uS, leak

% To eliminate currents, set them to zero:
%gNa = 0; %uS
%gCa = 0; % 10; %uS
%gK  = 0; % 70; %uS
%gNS = 0.05; % 0.5; %uS, non-specific. For some weird numerical

% reason, hard to set this to exactly zero.


R = 8314; %J 1000g/kg mol K
T = 295; %K

Na_o = 500; %mM
Na_i = 50; %mM
K_o = 10; %mM
Ca_o = 10; %mM

KMSI = 25e-6;  %mM, K_M of slow inward Ca current
KMCaP = 10^-3;  %mM, K_M of Ca pump
KMNS = 385e-6; %mM, K_M of non-specific current

%I_CaP = 15.0; %nA, for the calcium pump
I_CaP = 30;

ENa = 54; %mV, Nernst potential of Na
ECa = 65; %mV
EK = -75; %mV
ENS = -22; %mV
EL = 0; %mV

Cm = 17.5; %nF, membrane capacitance
Vol_i = 4.0; %nl, soma volume
F = 96500; %C/mol
kU = 100; %1/mM ms  % *** WATCH OUT *** kU AND kR ARE DEFINED ALSO IN OC_INF.m
kR = 0.238; %/ms    % 
B_i = 112.5e-3; %mM, internal concentration of cytosolic Ca buffer
KNaCa = 0.01;
DNaCa = 0.01;
gamma = 0.5;
r = 4;
nn = 4;  %DON'T NAME THIS THE SAME AS THE HH GATING PARAMETER!
INaK = 3.16; %nA, Na-K pump

