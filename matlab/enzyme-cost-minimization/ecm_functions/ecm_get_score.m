function [u_cost, u] = ecm_get_score(ecm_score,x,pp)

switch ecm_score,
  case 'mdf',      [u_cost, u] = ecm_mdf(x,pp);
  case 'memc2s',   [u_cost, u] = ecm_memc2s(x,pp);
  case 'memc2sp',  [u_cost, u] = ecm_memc2sp(x,pp);
  case 'emc1',     [u_cost, u] = ecm_emc1(x,pp);
  case 'emc2s',    [u_cost, u] = ecm_emc2s(x,pp);
  case 'emc2sp',   [u_cost, u] = ecm_emc2sp(x,pp);
  case 'emc3s',    [u_cost, u] = ecm_emc3s(x,pp);
  case 'emc3sp',   [u_cost, u] = ecm_emc3sp(x,pp);
  case 'emc4dm',   [u_cost, u] = ecm_emc4dm(x,pp);
  case 'emc4sm',   [u_cost, u] = ecm_emc4sm(x,pp);
  case 'emc4cm',   [u_cost, u] = ecm_emc4cm(x,pp);
  case 'emc4geom', [u_cost, u] = ecm_emc4geom(x,pp);
  otherwise,        u_cost = nan; u = nan * ones(nr,1);
end