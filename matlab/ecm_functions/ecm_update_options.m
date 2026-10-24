function ecm_options = ecm_update_options(network, ecm_options);

% ECM_UPDATE_OPTIONS - Helper function for ECM options (struct 'ecm_options')
%
% ecm_options = ecm_update_options(network, ecm_options)
%
% Update ECM options for a given network:
%  - adapt metabolite constraints
%  - insert protein cost weights

display('Updating the parameters for ECM');

% Insert default options, if not yet specified

ecm_options_default = ecm_default_options(network);
ecm_options         = join_struct(ecm_options_default,ecm_options);

% Adjust concentration constraints

conc_min = ecm_options.conc_min;
conc_max = ecm_options.conc_max;

ind_conc_min_equals_conc_max = find(conc_min==conc_max);
if length(ind_conc_min_equals_conc_max), 
  ecm_options.conc_fix(ind_conc_min_equals_conc_max) = ecm_options.conc_min(ind_conc_min_equals_conc_max);
  display('(ecm_update_options.m): metabolite bounds effectively constrain some metabolites to fixed values');
end

[nm,nr] = size(network.N);
if isempty(conc_min), conc_min = ecm_options.conc_min_default * ones(nm,1); end
if isempty(conc_max), conc_max = ecm_options.conc_max_default * ones(nm,1); end

% Replace cofactor concentrations

if length(ecm_options.replace_cofactors),  
  c_data = ecm_options.c_data;
  if find(isfinite(c_data)),
    display('  o Fixing some metabolite concentrations to data values');
    if 1,
      display('    (Using data value from the first state)');
      c_data_median = c_data(:,1);
    else
      display('    (Using median data value over states)');
      c_data_median = nanmedian(c_data,2);
    end    
    ll = label_names(ecm_options.replace_cofactors,network.metabolites);
    ll = ll(find(ll));
    ll = ll(find(isfinite(c_data_median(ll))));
    conc_min(ll) = c_data_median(ll);
    conc_max(ll) = c_data_median(ll);
  else
    display('  o No data values are given - cannot replace metabolite concentrations, as demanded, by data values');
  end
end

% Fix predefined concentrations

if length(ecm_options.fix_metabolites),
  display('  o Fixing some metabolite concentrations');
  for it = 1:length(ecm_options.fix_metabolites),
    ll = label_names(ecm_options.fix_metabolites(it),network.metabolite_names);
    conc_min(ll) = ecm_options.fix_metabolite_values(it);    
    conc_max(ll) = ecm_options.fix_metabolite_values(it);
  end
else
  display('  o All metabolite concentrations have finite ranges (no metabolites are fixed)');
end

ecm_options.conc_min = conc_min;
ecm_options.conc_max = conc_max;

if ~isfield(network,'metabolite_names'),
  network.metabolite_names = network.metabolites;
end
ind_conc_fix = find(ecm_options.conc_min == ecm_options.conc_max);
met_fix      = network.metabolite_names(ind_conc_fix);
conc_fix     = ecm_options.conc_min(ind_conc_fix);
if prod(size(met_fix))==0, met_fix={}; conc_fix=[]; end

ind_met_fix = label_names(met_fix, network.metabolite_names);

if setdiff(find(network.external),ind_met_fix),
  display('  o Some external metabolites in network will not be treated as fixed');
end

% 

ecm_options.conc_min    = conc_min;
ecm_options.conc_max    = conc_max;
ecm_options.conc_fix    = conc_fix;
ecm_options.met_fix     = met_fix;
ecm_options.ind_met_fix = ind_met_fix;

% -----------------------------
% adapt lambda_regularistion to typical magnitude of enzyme cost

median_cost = nanmedian(ecm_options.enzyme_cost_weights);

if isnan(median_cost), median_cost = 1; end 

if length(ecm_options.lambda_reg_factor),
  ecm_options.lambda_regularisation = ecm_options.lambda_reg_factor * median_cost;
end

% Adapt print_graphics

if ecm_options.show_graphics == 0, ecm_options.print_graphics = 0;         end
if ecm_options.compute_tolerance,  ecm_options.tolerance_from_hessian = 0; end

% -----------------------------

if ~isfield(ecm_options,'N_subunits_in_complex'),
  ecm_options.N_subunits_in_complex  = ones(size(ecm_options.enzyme_cost_weights));
end

if ~isfield(ecm_options,'N_cat_sites_in_complex'),
  ecm_options.N_cat_sites_in_complex = ones(size(ecm_options.enzyme_cost_weights));
end

ecm_options.ecm_scores = column(ecm_options.ecm_scores);