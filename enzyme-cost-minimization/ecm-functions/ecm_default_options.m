function ecm_options = ecm_default_options(network, model_name)

% ECM_DEFAULT_OPTIONS -  Default options for ECM algorithm
%
% ecm_options = ecm_default_options(network, model_name)
%
% The default values are read from the SBtab Config file 'enzyme-cost-minimization_Config.tsv'
% (in the main code folder 'enzyme-cost-minimization'; TableID 'ConfigureECM'). To change the
% defaults, edit this file. The ECM options also contain all options of parameter balancing
% (see 'parameter_balancing_options'); options given in both files take their values from the ECM file.
% Options that depend on the network (ind_scored_enzymes, enzyme_cost_weights, show_metabolites)
% are set by this function if they are empty in the Config file.
%
% ECM options:
%
% General options and model name
% ------------------------------
%   ecm_options.verbose                (Boolean)
%   ecm_options.model_name             (string) model_name (function argument; if not given, taken from the Config file)
%   ecm_options.model_id               (string) model id (for file names)
%   ecm_options.run_id                 (string) simulation run id (for file names)
%
% Metabolite constraints
% ----------------------
%   ecm_options.conc_min               (float vector)      minimal metabolite concentrations (mM), with nans
%   ecm_options.conc_max               (float vector)      maximal metabolite concentrations (mM), with nans
%   ecm_options.conc_min_default       (float)             default minimal metabolite concentration (mM)
%   ecm_options.conc_max_default       (float)             default maximal metabolite concentration (mM)
%   ecm_options.fix_metabolites        (Boolean vector)    metabolites with fixed concentrations (overrides values from input model)
%   ecm_options.fix_metabolite_values  (vector)            fixed metabolite concentrations (see ecm_options.fix_metabolites)
%   ecm_options.met_fix                                    metabolites with fixed concentrations
%   ecm_options.ind_met_fix                                indices of metabolites with fixed concentrations
%   ecm_options.conc_fix               (float vector)      fixed metabolite concentrations (mM) (see ecm_options.met_fix)
%   ecm_options.replace_cofactors      (string cell array) names of cofactors; this is only used in ecm_update_options.m -
%                                                          bounds for cofactor concentrations are automatically set based on
%                                                          median value from concentration data
% Kinetic data
% ------------
%   ecm_options.kinetic_data                     matlab struct with kinetic data; see data structures in mnt toolbox
%   ecm_options.reaction_column_names            reaction column names (in data file) for loading of kinetic data
%   ecm_options.compound_column_names            compound column names (in data file) for loading of kinetic data
%   ecm_options.KM_lower               (float)   lower bound on KM values (mM)
%   ecm_options.Keq_upper              (float)   upper bound on Keq values
%   ecm_options.flag_given_kinetics    (Boolean) use given kinetics
%   ecm_options.kcat_usage             (string)  use kcat data values ('use')
%   ecm_options.kcat_upper             (float)   upper bound on Kcat values (1/s)
%   ecm_options.kcat_lower             (float)   lower bound on Kcat values (1/s)
%   ecm_options.kcatr_lower            (float)   lower bound on reverse Kcat values (1/s)
%   ecm_options.kcat_prior_median      (float)   median for Kcat values (1/s)
%   ecm_options.kcat_prior_log10_std   (float)   std deviation of log10 Kcat values
%
% Parameter balancing
% -------------------
%   ecm_options.use_pseudo_values      (Boolean) use pseudo values during parameter balancing
%   ecm_options.GFE_fixed              (Boolean) keep GFE values fixed during parameter balancing
%   ecm_options.insert_Keq_from_data   (Boolean) insert again Keq values from data after parameter balancing
%   ecm_options.n_samples              (integer) number of samples in posterior sampling
%
% Enzyme cost weights
% -------------------
%   ecm_options.ind_scored_enzymes     (integer vector) indices of reactions with enzymes to be scored in cost function (empty: all reactions)
%   ecm_options.enzyme_cost_weights    (float vector)   numerical enzyme cost weights (empty: ones)
%   ecm_options.use_cost_weights       (string)         usage of cost weights
%   ecm_options.parameter_prior_file   (string)         parameter prior file filename, optional
%
% ECM calculation
% ---------------
%   ecm_options.ignore_regulation           (Boolean)      disregard enzyme regulation
%   ecm_options.initial_choice              (string)       method for computing inital solution {'mdf' 'polytope_center', 'interval_center',
%                                                          'given_x_start'}
%   ecm_options.x_start                     (float vector) initial solution vector, only for ecm_options.initial_choice='given_x_start'
%   ecm_options.multiple_conditions_n       (int)          number of experimental conditions to be considered
%   ecm_options.multiple_conditions_anticipate  (boolean)  optimise a single enzyme profile for multiple conditions?
%   ecm_options.multiple_starting_points    (boolean)      use multiple starting points
%   ecm_options.ecm_scores                  (string cell array) list of ecm scores to be used; see 'help ecm_scores'
%   ecm_options.lambda_regularisation       (float)        coefficient for metabolite concentration regularisation term
%   ecm_options.lambda_reg_factor           (float)        regularisation factor for metabolite concentration regularisation term
%                                                          (parameter for choosing ecm_options.lambda_regularisation automatically, in ecm_update_options; set to [] in order to not use it)
%   ecm_options.use_linear_cost_constraints (Boolean)      tighter constraints derived upper bound on enzyme cost
%   ecm_options.fluctuations_safety_margin  (float)        safety margin (# std dev) to counter protein number fluctuations
%   ecm_options.cell_volume                 (float)        cell volume in m^3 (needed for safety margin; default value for E coli)
%   ecm_options.maximal_u_cost              (float)        used with option "use_linear_cost_constraints"
%   ecm_options.compute_hessian             (Boolean)      compute Hessian matrix
%   ecm_options.compute_elasticities        (Boolean)      compute elasticities
%   ecm_options.compute_tolerance           (Boolean)      compute tolerance ranges
%   ecm_options.cost_tolerance_factor       (float)        factor defining the cost tolerance
%   ecm_options.tolerance_from_hessian      (Boolean)      compute tolerances from the Hessian matrix
%   ecm_options.fix_Keq_in_sampling         (Boolean)      keep Keq values fixed in sampling
%   ecm_options.fix_thermodynamics_by_adjusting_Keq  (Boolean) enforce thermodynamically feasibile fluxes by adjusting Keq
%
% Metabolic data (for comparison)
% -------------------------------
%   ecm_options.c_data  (float vector) metabolite concentration data, with nans
%   ecm_options.u_data  (float vector) enzyme concentration data, with nans
%
% Graphics
% --------
%   ecm_options.print_graphics         (Boolean)           save graphics to files
%   ecm_options.show_graphics          (Boolean)           show graphics
%   ecm_options.show_metabolites       (string cell array) metabolites shown in graphics (empty: network.metabolites)
%   ecm_options.metabolite_order_file  (string)            filename of file containing metabolite order for graphics
%   ecm_options.reaction_order_file    (string)            filename of file containing reaction order for graphics

% Argument 'model_name' is optional

eval(default('network','[]'));

if isempty(network),
  network = network_construct;%(N,reversible,ind_external,metabolites,actions,no_graphics,regulation_matrix,flag_kinetics)
end

% parameter balancing options, overridden by the ECM options from the Config file

config_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'enzyme-cost-minimization_Config.tsv');
ecm_options = join_struct(parameter_balancing_options, load_config_table(config_file));

if exist('model_name','var'),
  ecm_options.model_name = model_name;
end

% options that depend on the network

if isempty(ecm_options.ind_scored_enzymes),
  ecm_options.ind_scored_enzymes  = [1:length(network.actions)]';
end
if isempty(ecm_options.enzyme_cost_weights),
  ecm_options.enzyme_cost_weights = ones(length(ecm_options.ind_scored_enzymes),1);
end
if isempty(ecm_options.show_metabolites),
  ecm_options.show_metabolites    = network.metabolites;
end
