function [u_tot, u] = measures_for_enzyme_costs_score_fsc4cmr(x,v,network,ind_scored_enzymes,enzyme_cost_weights)

% [u_tot, u] = measures_for_enzyme_costs_score_fsc4cmr(x,v,network,ind_scored_enzymes,enzyme_cost_weights)

delta_G_by_RT =  sign(v) .* [network.N' * x - log(network.kinetics.Keq)];

network.kinetics.type = 'ms';
network.kinetics.u    = ones(size(network.kinetics.u)); 

w = network_velocities(exp(x),network);

u = v./w;

u_tot = sum(enzyme_cost_weights.*u(ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end
