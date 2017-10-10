function [u_cost, u, w] = ecm_fitting_3sp(x,pp)

% [u_cost, u] = ecm_fitting(x,pp)
% sum of squared residuals for data fit (metabolite and enzyme levels, both on log scale), assuming emc3sp formula for enzyme demand

delta_G_by_RT = pp.N_forward' * x - pp.log_Keq_forward;

[~, u, w] = ecm_emc3sp(x,pp);

if sum(delta_G_by_RT>0),
  u_cost = 10^20 * sum(delta_G_by_RT(find(delta_G_by_RT>0)));
end

ind_c = find(isfinite(pp.ln_c_data));
ind_u = find(isfinite(pp.u_data));

u_cost = sum([[x(ind_c)-pp.ln_c_data(ind_c)] ./ pp.ln_c_std(ind_c)].^2) ...
       + sum([[u(ind_u)-pp.u_data(ind_u)]    ./ pp.u_std(ind_u)   ].^2)
