function [u_cost, u] = ecm_get_score(fsc_score,x,pp)

switch fsc_score,
  case 'mdf',      [u_cost, u] = ecm_mdf(x,pp);
  case 'mfsc2sub', [u_cost, u] = ecm_min_ecf2s(x,pp);
  case 'mfsc2',    [u_cost, u] = ecm_min_ecf2sp(x,pp);
  case 'ecf1',     [u_cost, u] = ecm_ecf1(x,pp);
  case 'ecf2s',    [u_cost, u] = ecm_ecf2s(x,pp);
  case 'ecf2sp',   [u_cost, u] = ecm_ecf2sp(x,pp);
  case 'ecf3s',    [u_cost, u] = ecm_ecf3s(x,pp);
  case 'ecf3sp',   [u_cost, u] = ecm_ecf3sp(x,pp);
  case 'ecf4dmr',  [u_cost, u] = ecm_ecf4dmr(x,pp);
  case 'ecf4smr',  [u_cost, u] = ecm_ecf4smr(x,pp);
  case 'ecf4cmr',  [u_cost, u] = ecm_ecf4cmr(x,pp);
  case 'ecf4geom', [u_cost, u] = ecm_ecf4geom(x,pp);
  otherwise,        u_cost = nan; u = nan * ones(nr,1);
end