function [u_cost, u, w] = ecm_emc2s(x,pp)

% [u_cost, u] = ecm_emc2s(x,pp)

delta_G_by_RT = pp.N_forward' * x - pp.log_Keq_forward;

w = [ pp.kc_forward .* [1 - exp(delta_G_by_RT)] ];
u = abs(pp.v) ./ w;
u(pp.v==0) = 0;

if sum(delta_G_by_RT>0),
  u(find(delta_G_by_RT>0)) = 10^20 * max(delta_G_by_RT);
  warning('Unfeasible metabolite profile');
end

u_cost = sum(pp.enzyme_cost_weights .* u(pp.ind_scored_enzymes));
