function rates = ecm_get_specific_rates(fsc_score,x,pp,ind)

switch fsc_score,
  case 'mdf',      [u_tot, u, rates] = ecm_mdf(x,pp);
  case 'mfsc2sub', [u_tot, u, rates] = ecm_min_ecf2s(x,pp);
  case 'mfsc2',    [u_tot, u, rates] = ecm_min_ecf2sp(x,pp);
  case 'ecf1',     [u_tot, u, rates] = ecm_ecf1(x,pp);
  case 'ecf2s',    [u_tot, u, rates] = ecm_ecf2s(x,pp);
  case 'ecf2sp',   [u_tot, u, rates] = ecm_ecf2sp(x,pp);
  case 'ecf3s',    [u_tot, u, rates] = ecm_ecf3s(x,pp);
  case 'ecf3sp',   [u_tot, u, rates] = ecm_ecf3sp(x,pp);
  case 'ecf4dmr',  [u_tot, u, rates] = ecm_ecf4dmr(x,pp);
  case 'ecf4smr',  [u_tot, u, rates] = ecm_ecf4smr(x,pp);
  case 'ecf4cmr',  [u_tot, u, rates] = ecm_ecf4cmr(x,pp);
  case 'ecf4geom', [u_tot, u, rates] = ecm_ecf4geom(x,pp);
  otherwise,       u_tot = nan; u = nan * ones(nr,1); rates = nan * ones(nr,1);
end

if exist('ind','var'),  rates = rates(ind);  end