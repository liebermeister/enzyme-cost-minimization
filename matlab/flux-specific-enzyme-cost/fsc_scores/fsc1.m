function [u_tot, u] = measures_for_enzyme_costs_score_fsc1(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,ind_scored_enzymes,enzyme_cost_weights)

% [u_tot, u] = measures_for_enzyme_costs_score_fsc1(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,ind_scored_enzymes,enzyme_cost_weights)

u = abs(v) ./ kc_forward;

u_tot = sum(enzyme_cost_weights .* u(ind_scored_enzymes));