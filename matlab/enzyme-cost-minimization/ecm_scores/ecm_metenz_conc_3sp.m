function [u_cost, u, w] = ecm_metenz_conc_3sp(x,pp)

% [u_cost, u] = ecm_metenz_conc_3sp(x,pp)
% sum of total metabolite and enzyme levels, assuming emc3sp formula for enzyme demand

delta_G_by_RT = pp.N_forward' * x - pp.log_Keq_forward;

[~, u, w] = ecm_emc3sp(x,pp);

if sum(delta_G_by_RT>0),
  u_cost = 10^20 * sum(delta_G_by_RT(find(delta_G_by_RT>0)));
end

u_cost = sum(exp(x)) + sum(u);
