function [u_tot, u, w] = ecm_emc1(x,pp)

% [u_tot, u] = ecm_emc1(x,pp)

u = abs(pp.v) ./ pp.kc_forward;

u_tot = sum(pp.enzyme_cost_weights .* u(pp.ind_scored_enzymes));

w = pp.kc_forward;
