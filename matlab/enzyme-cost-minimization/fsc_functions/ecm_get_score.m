function [my_u_cost, my_u] = ecm_get_score(fsc_score,my_x,pp)

switch fsc_score,
  case 'mdf',      [my_u_cost, my_u] = ecm_mdf(my_x,pp);
  case 'mfsc2sub', [my_u_cost, my_u] = ecm_min_ecf2s(my_x,pp);
  case 'mfsc2',    [my_u_cost, my_u] = ecm_min_ecf2sp(my_x,pp);
  case 'ecf1',     [my_u_cost, my_u] = ecm_ecf1(my_x,pp);
  case 'ecf2s',    [my_u_cost, my_u] = ecm_ecf2s(my_x,pp);
  case 'ecf2sp',   [my_u_cost, my_u] = ecm_ecf2sp(my_x,pp);
  case 'ecf3s',    [my_u_cost, my_u] = ecm_ecf3s(my_x,pp);
  case 'ecf3sp',   [my_u_cost, my_u] = ecm_ecf3sp(my_x,pp);
  case 'ecf4dmr',  [my_u_cost, my_u] = ecm_ecf4dmr(my_x,pp);
  case 'ecf4smr',  [my_u_cost, my_u] = ecm_ecf4smr(my_x,pp);
  case 'ecf4cmr',  [my_u_cost, my_u] = ecm_ecf4cmr(my_x,pp);
  case 'ecf4geom', [my_u_cost, my_u] = ecm_ecf4geom(my_x,pp);
  otherwise,        my_u_cost = nan; my_u = nan * ones(nr,1);
end