clear all; close all; clc;
all_constants;

T = 20;
Istim=0;
V0=0;
Ca=0;

dt = 1e-5;
%time = 0:dt:T;
time = T;

%%
[time, V, Ca_i, O_c, INa, ICa, IISI, INS, IL, IK, IR, ICaP, INaCa] = Semireduced15sim(time, Istim, V0, Ca);


%%
dV =  - (INa + ICa + IISI + INS + IK + IR + IL + INaCa + INaK + ICaP ) / Cm;
figure()
plot(V, dV)

%% Save output
save('../burst_data.mat', 'V', 'dV', 'INa', 'ICa', 'IISI', 'INS', 'IK', 'IR', ...
    'IL', 'INaCa', 'INaK', 'ICaP', 'Cm', 'time')