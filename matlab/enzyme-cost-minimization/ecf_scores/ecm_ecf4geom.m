function [u_tot, u] = ecm_ecf4geom(x,pp)

% [u_tot, u] = ecm_ecf4geom(x,pp)
% geometric mean of cs and ds rate law

delta_G_by_RT =  pp.N_forward' * x - pp.log_Keq_forward;

[u_tot, u_ds] = ecm_ecf3sp(x,pp);
[u_tot, u_cs] = ecm_ecf4cmr(x,pp);

u = sqrt(u_ds .* u_cs);

u_tot = sum(pp.enzyme_cost_weights .* u(pp.ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end
