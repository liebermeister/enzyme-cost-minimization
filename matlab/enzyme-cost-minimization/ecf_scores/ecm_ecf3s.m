function [u_tot, u] = ecm_ecs3s(x,pp)

% [u_tot, u] = ecm_ecs3s(x,pp)

delta_G_by_RT =  pp.N_forward' * x - pp.log_Keq_forward;

u = abs(pp.v) ./ [ pp.kc_forward .* [1 - exp(delta_G_by_RT)] ./ [1 + pp.km_forward .* exp(-pp.M_forward * x)]  ];

u_tot = sum(pp.enzyme_cost_weights .* u(pp.ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end
