function dvdt = v_cline(V, Ca_i, Istim)

dYdt = reduced15neuron(10, [V;Ca_i], Istim);

dvdt = dYdt(1);

