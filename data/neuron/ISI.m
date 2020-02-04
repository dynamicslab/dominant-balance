function ii = ISI(s, V, Ca_i)

all_constants;

% ii = gSI * s .* (V - ECa) .* (0.11/600)*KMSI./(KMSI.^2 + Ca_i.^2);
ii = gSI * s .* (V - ECa) .* KMSI./(KMSI + Ca_i);
