function [u_tot, u, w] = ecm_ecf3s(x,pp)

% [u_tot, u] = ecm_ecf3s(x,pp)

delta_G_by_RT =  pp.N_forward' * x - pp.log_Keq_forward;

w = [ pp.kc_forward .* [1 - exp(delta_G_by_RT)] ./ [1 + pp.km_forward .* exp(-pp.M_forward * x)]  ];
u = abs(pp.v) ./ w;

u_tot = sum(pp.enzyme_cost_weights .* u(pp.ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end
