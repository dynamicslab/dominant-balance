function [time, V, Ca_i, O_c, INa, ICa, IISI, INS, IL, IK, IR, ICaP, INaCa] = Semireduced15sim(time, Istim, V, Ca)

%INITIAL CONDITIONS
V; %initial membrane potential (mV)
%Ca_i = 0;
Ca_i = Ca*10^-6; %mM


%INTEGRATE
options = odeset('RelTol',1e-4,'AbsTol',1e-5);
%[T,Y] = ode15s(@(t, y) semireduced15neuron(t, y, Istim),[0 time*1000], ...
%  [V ; Ca_i ; 0.15; minf(V) ; hinf(V) ; ninf(V); linf(V)],options);
[T,Y] = ode15s(@(t, y) semireduced15neuron(t, y, Istim),[0, time*1000], ...
  [V ; Ca_i ; 0.15; minf(V) ; hinf(V) ; ninf(V); linf(V)],options);

V    = Y(:,1);
Ca_i = Y(:,2);
time = T/1000;
O_c   = Y(:,3);
m    = Y(:,4);
h    = Y(:,5);
d    = dinf(V);
f    = finf(V);
s    = sinf(V);
b    = binf(V);
n    = Y(:,6);
l    = Y(:,7);


figure(6)
subplot(3,2,[1 2])
plot(T/1000, V)
xlabel('Time (sec)');
ylabel('Membrane Potential (mV)');
subplot(3,2,3)
plot(T, Ca_i,'k')
legend('[Ca]_i');
ylim([0 600e-6]);
subplot(3,2,4)
plot(T/1000, O_c,'k')
legend('O_c');
ylim([0 1]);
subplot(3,2,[5 6])
plot(T/1000, [m h d f s b n l])
legend('m', 'h', 'd', 'f', 's', 'b', 'n', 'l');
ylim([0 1]);

gamma = 0; all_constants;
INa = gNa * m.^3 .* h .* (V - ENa);
ICa = gCa * d.^2 .* f .* (Y(:,1) - ECa);
% ISI = gSI * Y(:,8) .* (V - ECa) .* (0.11/600)*KMSI./(KMSI.^2 + Ca_i.^2);
IISI = ISI(s, V, Ca_i);
INS = gNS * b .* (V - ENS) .* Ca_i./(Ca_i + KMNS);
IL = gL * (V - EL);
IK = gK * n.^4 .* l .* (V - EK);
Z = 1;  %Z was not defined in the paper; I assume it's a charge.
IR = gR * (V - EK + 5.66) ./ (1 + exp((V - EK - 15.3)*Z * F/(R*T)));



ICaP = I_CaP * (Ca_i ./ (Ca_i + KMCaP));
DFin = (Na_i)^r * (Ca_o) * exp((r-2) * gamma * V *F/ (R*T));
DFout = (Na_o)^r * (Ca_i) .* exp((r-2) * (gamma-1) * V * F/ (R*T));
S = 1 + DNaCa * (Ca_i.*(Na_o)^r + Ca_o*(Na_i)^r);
INaCa = KNaCa * (DFin - DFout) ./ S;