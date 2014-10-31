function fsc_options_default = fsc_default_options(network)

% fsc_options_default = fsc_default_options(network)

[nm,nr] = size(network.N);

fsc_options_default.model_name               = 'model';
fsc_options_default.run_id                   = 'run';
fsc_options_default.network_CoHid            = network;

% list of metabolite names with fixed concentrations (overrides values from input model)
fsc_options_default.fix_metabolites          = {}; 
fsc_options_default.fix_metabolite_values    = [];

fsc_options_default.insert_original_equilibrium_constants = 0;

fsc_options_default.kcat_usage               = 'use';
fsc_options_default.kcat_prior_median        = 350;   % similar to median in glycolysis+tca
fsc_options_default.kcat_prior_log10_std     = 0.1;   % reduce spread of kcat values
fsc_options_default.kcat_lower               = 50;   % 
fsc_options_default.kcat_upper               = 2000;  % 

% promlematic if some concentrations are thermodynamically forced to be very low
fsc_options_default.KM_lower                 = [];  % mM
fsc_options_default.Keq_upper                = [];

fsc_options_default.conc_min_default         = 10^-3; % 0.001 mM
fsc_options_default.conc_max_default         = 10;  % 10 mM
fsc_options_default.conc_min                 = nan*ones(nm,1);
fsc_options_default.conc_max                 = nan*ones(nm,1);

fsc_options_default.c_data                   = [];
fsc_options_default.u_data                   = [];
fsc_options_default.kinetic_data             = [];

% column names (in data file) for loading of kinetic data
fsc_options_default.reaction_column_names = [];
fsc_options_default.compound_column_names = [];

fsc_options_default.flag_given_kinetics      = 0;
fsc_options_default.enzyme_cost_weights      = [];

fsc_options_default.lambda_regularisation    =  10^-5; 
fsc_options_default.fsc_scores               = {'ecf1'};
fsc_options_default.initial_choice           = 'mdf'; 
fsc_options_default.quantity_info_file       = [];
fsc_options_default.ind_scored_enzymes       = 1:length(network.actions);

fsc_options_default.variation_for_relaxed_optimality = 1;
fsc_options_default.multiple_starting_points = 0;

fsc_options_default.print_graphics           = 0;
fsc_options_default.show_graphics            = 1;
fsc_options_default.show_metabolites         = network.metabolites;
