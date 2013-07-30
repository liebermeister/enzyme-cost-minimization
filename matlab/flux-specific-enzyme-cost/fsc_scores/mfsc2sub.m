function [f, u] = measures_for_enzyme_costs_score_mfsc2sub(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)

% [f, u] = measures_for_enzyme_costs_score_mfsc2sub(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)

delta_G_by_RT = N_forward' * x - log_Keq_forward;

u = abs(v) ./ [ kc_forward .* [1 - exp(delta_G_by_RT)] ];

f = max(enzyme_cost_weights .* u(ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  f = 10^20 * max(delta_G_by_RT);
end
