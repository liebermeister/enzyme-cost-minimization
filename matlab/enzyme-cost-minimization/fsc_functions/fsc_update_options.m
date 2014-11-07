function fsc_options = fsc_update_options(network, fsc_options, r_orig);

% fsc_options = fsc_update_options(network, fsc_options, r_orig)
%
% Update ECM options for a given network:
%  - adapt metabolite constraints
%  - insert protein cost weights
%  - adjust Keq values


fsc_options_default = fsc_default_options(network);
fsc_options         = join_struct(fsc_options_default,fsc_options);

% Adjust concentration constraints

conc_min = fsc_options.conc_min;
conc_max = fsc_options.conc_max;
c_data   = fsc_options.c_data  ;

[nm,nr] = size(network.N);
if isempty(conc_min), conc_min = fsc_options.conc_min_default * ones(nm,1); end
if isempty(conc_max), conc_max = fsc_options.conc_max_default * ones(nm,1); end

% Replace cofactor concentrations

display('  Replacing cofactor concentrations');

c_data_median = nanmedian(c_data,2);
ll = label_names(fsc_options.replace_cofactors,network.metabolites);
ll = ll(find(ll));
ll = ll(find(isfinite(c_data_median(ll))));
conc_min(ll) = c_data_median(ll);
conc_max(ll) = c_data_median(ll);

% Fix predefined concentrations

display('  (possibly) fixing some metabolite concentrations');
for it = 1:length(fsc_options.fix_metabolites),
  ll = label_names(fsc_options.fix_metabolites(it),network.metabolites);
  conc_min(ll) = fsc_options.fix_metabolite_values(it);
  conc_max(ll) = fsc_options.fix_metabolite_values(it);
end

ind_conc_fix = find(fsc_options.conc_min == fsc_options.conc_max);
met_fix      = network.metabolites(ind_conc_fix);
conc_fix     = fsc_options.conc_min(ind_conc_fix);

ind_met_fix = label_names(met_fix, network.metabolites);

if setdiff(find(network.external),ind_met_fix),
  display('  Some external metabolites in network will not be treated as fixed');
end

% 

fsc_options.conc_min    = conc_min;
fsc_options.conc_max    = conc_max;
fsc_options.conc_fix    = conc_fix;
fsc_options.met_fix     = met_fix;
fsc_options.ind_met_fix = ind_met_fix;


% -------------------------------------------
% set enzyme cost weights

switch fsc_options.use_cost_weights,
  case 'none',
    fsc_options.enzyme_cost_weights = ones(length(fsc_options.ind_scored_enzymes),1);
  case 'protein_size',
    fsc_options.enzyme_cost_weights = network.enzyme_size(fsc_options.ind_scored_enzymes);
  case 'aa_composition',
    fsc_options.enzyme_cost_weights = network.akashi_protein_cost(fsc_options.ind_scored_enzymes);
end  

if length(fsc_options.enzyme_cost_weights),
  enzyme_cost_weights= fsc_options.enzyme_cost_weights / median(fsc_options.enzyme_cost_weights);
else,
  enzyme_cost_weights= ones(length(ind_scored_enzymes),1);
end


% -----------------------------
% if desired, insert original Keq values

if fsc_options.insert_Keq_from_data,
  display('    Using predefined equilibrium constants exactly'); 
  r.Keq = r_orig.Keq;
end


% Adapt print_graphics

if fsc_options.show_graphics == 0, fsc_options.print_graphics = 0;         end
if fsc_options.compute_tolerance,  fsc_options.tolerance_from_hessian = 0; end
