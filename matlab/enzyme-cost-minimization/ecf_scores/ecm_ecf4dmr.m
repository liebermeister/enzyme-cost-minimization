function [u_tot, u, w] = ecm_ecf4cmr(x,pp)

% [u_tot, u] = ecm_ecf4cmr(x,pp)

network = pp.network;

delta_G_by_RT =  sign(pp.v) .* [network.N' * x - log(pp.network.kinetics.Keq)];

network.kinetics.type = 'ds';
network.kinetics.u    = ones(size(network.N,2),1); 

w = network_velocities(exp(x),network);

u = pp.v./w;

u_tot = sum(pp.enzyme_cost_weights.* u(pp.ind_scored_enzymes));

if sum(delta_G_by_RT>0),
  u_tot = 10^20*max(delta_G_by_RT);
end

% sprintf('%f',u_tot)
