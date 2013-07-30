function [f,u] = fsc_obdw(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)

% [f,u] = fsc_mtdf(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)

delta_G_by_RT = N_forward' * x - log_Keq_forward;

f = max( abs(1./v) .* delta_G_by_RT(ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  f = 10^20*max(delta_G_by_RT);
end

u = nan * v;
