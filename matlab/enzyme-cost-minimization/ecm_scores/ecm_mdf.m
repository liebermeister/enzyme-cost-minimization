function [f,u,w] = ecm_mdf(x,pp)

% [f,u] = ecm_mdf(x,pp)

delta_G_by_RT = pp.N_forward' * x - pp.log_Keq_forward;

f = max(delta_G_by_RT(pp.ind_scored_enzymes));

if sum(delta_G_by_RT > 0),
  f = 10^20*max(delta_G_by_RT);
end

u = nan * pp.v;
w = nan * pp.v;
