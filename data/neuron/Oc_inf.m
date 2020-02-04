function Oc = Oc_inf(Ca_i)

kU = 100; %1/mM ms
kR = 0.238; %/ms

Oc = kU*Ca_i / (kR + kU*Ca_i);


