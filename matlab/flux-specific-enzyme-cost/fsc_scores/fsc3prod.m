function [u_tot, u] = measures_for_enzyme_costs_score_fsc3prod(x,v,M_forward,Mprod_forward,N_forward,log_Keq_forward,kc_forward,km_forward,kmprod_forward,ind_scored_enzymes,enzyme_cost_weights)

% [u_tot, u] = measures_for_enzyme_costs_score_fsc3prod(x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)

delta_G_by_RT =  N_forward' * x - log_Keq_forward;

w =  [ kc_forward .* [1 - exp(delta_G_by_RT)] ./ [[  km_forward  ./ exp(abs(M_forward) * x)] .* [1 + 1./km_forward .* exp(abs(M_forward) * x) + 1./kmprod_forward .* exp(abs(Mprod_forward) * x)]  ]];

u = abs(v) ./ w;

u_tot = sum(enzyme_cost_weights .* u(ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end

%sprintf('%f',u_tot)

