% ECM_GET_SCORE - Helper function
% 
% function [u_cost, u_preemptive] = ecm_get_score(ecm_score,x,pp)

function [u_cost, u, u_preemptive] = ecm_get_score(ecm_score,x,pp)

switch ecm_score,
  case 'mdf',      [u_cost, u] = ecm_mdf(x,pp);      % only max Delta G counts
  case 'memc2s',   [u_cost, u] = ecm_memc2s(x,pp);   % max values
  case 'memc2sp',  [u_cost, u] = ecm_memc2sp(x,pp);
  case 'emc1',     [u_cost, u] = ecm_emc1(x,pp);     % sum values
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

% Support for multiple conditions
% Assume that all reaction vectors are actually a concatenation of several (condition-specific) vectors

if pp.multiple_conditions,
  display(sprintf('  Running multi-condition optimisation with %d conditions',pp.multiple_conditions_n));
  switch ecm_score,
    case {'emc1', 'emc2s', 'emc2sp', 'emc3s',  'emc3sp', 'emc4dm', 'emc4sm', 'emc4cm', 'emc4geom'},
      %% Recompute enzyme values and total cost
      nu                 = length(u) / pp.multiple_conditions_n;
      u_preemptive       = nanmax(reshape(u,nu,pp.multiple_conditions_n),[],2);
      ind_scored_enzymes = pp.ind_scored_enzymes(find(pp.ind_scored_enzymes <= nu));
      u_cost             = sum( pp.enzyme_cost_weights(ind_scored_enzymes) .* u_preemptive(ind_scored_enzymes) );
    otherwise, 
      error('With this cost score, multiple conditions are not supported '); 
  end
else
  u_preemptive = [];
end