% ECM_GET_SCORE - Helper function
% 
% function [u_cost, u_preemptive] = ecm_get_score(ecm_score,x,pp)
%
% Possible cost functions: 
%  Max-min driving force
%   mdf
%   
%  Enzyme cost based on different rate law approximations (see Noor et al 2016)
%   emc1
%   emc2s
%   emc2sp        
%   emc3s         
%   emc3sp        
%   emc4dm        
%   emc4sm        
%   emc4cm        
%   emc4geom      
%   
%  Maximal cost among enzymes
%   max-1
%   max-2s        
%   max-2sp
%   max-3s        
%   max-3sp        
%   max-4cm 
%   
%  Optimal fit (sum of squared residuals) for metabolite and  enzymes data
%   fitting-1     
%   fitting-2s    
%   fitting-2sp   
%   fitting-3s    
%   fitting-3sp   
%   fitting-4cm   
%   
%  Sum of metabolite and enzyme concentrations
%   metenzconc-1  
%   metenzconc-2s 
%   metenzconc-2sp
%   metenzconc-3s 
%   metenzconc-3sp
%   metenzconc-4cm
%   
%  Sum of mass-weighted metabolite and enzyme concentrations
%   metenzmass-1  
%   metenzmass-2s 
%   metenzmass-2sp
%   metenzmass-3s 
%   metenzmass-3sp
%   metenzmass-4cm

function [u_cost, u, u_preemptive] = ecm_get_score(ecm_score,x,pp)

switch ecm_score,
  case 'mdf',            [u_cost, u] = ecm_mdf(x,pp);      % only max Delta G counts  
  case 'emc1',           [u_cost, u] = ecm_emc1(x,pp);     % sum values
  case 'emc2s',          [u_cost, u] = ecm_emc2s(x,pp);
  case 'emc2sp',         [u_cost, u] = ecm_emc2sp(x,pp);
  case 'emc3s',          [u_cost, u] = ecm_emc3s(x,pp);
  case 'emc3sp',         [u_cost, u] = ecm_emc3sp(x,pp);
  case 'emc4dm',         [u_cost, u] = ecm_emc4dm(x,pp);
  case 'emc4sm',         [u_cost, u] = ecm_emc4sm(x,pp);
  case 'emc4cm',         [u_cost, u] = ecm_emc4cm(x,pp);
  case 'emc4geom',       [u_cost, u] = ecm_emc4geom(x,pp);  
  case 'max-1',          [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc1');   % max values
  case 'max-2s',         [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc2s');
  case 'max-2sp',        [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc2sp');
  case 'max-3s',         [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc3s');
  case 'max-3sp',        [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc3sp');
  case 'max-4cm',        [u_cost, u] = ecm_max_enzyme_cost(x,pp,'emc4cm');
  case 'fitting-1',      [u_cost, u] = ecm_fitting(x,pp,'emc1');
  case 'fitting-2s',     [u_cost, u] = ecm_fitting(x,pp,'emc2s');
  case 'fitting-2sp',    [u_cost, u] = ecm_fitting(x,pp,'emc2sp');
  case 'fitting-3s',     [u_cost, u] = ecm_fitting(x,pp,'emc3s');
  case 'fitting-3sp',    [u_cost, u] = ecm_fitting(x,pp,'emc3sp');
  case 'fitting-4cm',    [u_cost, u] = ecm_fitting(x,pp,'emc4cm');
  case 'metenzconc-1',   [u_cost, u] = ecm_metenz_conc(x,pp,'emc1');
  case 'metenzconc-2s',  [u_cost, u] = ecm_metenz_conc(x,pp,'emc2s');
  case 'metenzconc-2sp', [u_cost, u] = ecm_metenz_conc(x,pp,'emc2sp');
  case 'metenzconc-3s',  [u_cost, u] = ecm_metenz_conc(x,pp,'emc3s');
  case 'metenzconc-3sp', [u_cost, u] = ecm_metenz_conc(x,pp,'emc3sp');
  case 'metenzconc-4cm', [u_cost, u] = ecm_metenz_conc(x,pp,'emc4cm');
  case 'metenzmass-1',   [u_cost, u] = ecm_metenz_mass(x,pp,'emc1');
  case 'metenzmass-2s',  [u_cost, u] = ecm_metenz_mass(x,pp,'emc2s');
  case 'metenzmass-2sp', [u_cost, u] = ecm_metenz_mass(x,pp,'emc2sp');
  case 'metenzmass-3s',  [u_cost, u] = ecm_metenz_mass(x,pp,'emc3s');
  case 'metenzmass-3sp', [u_cost, u] = ecm_metenz_mass(x,pp,'emc3sp');
  case 'metenzmass-4cm', [u_cost, u] = ecm_metenz_mass(x,pp,'emc4cm');
  
  otherwise,             error('Unknown ecm_score');
end

u(pp.v==0) = 0;

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
