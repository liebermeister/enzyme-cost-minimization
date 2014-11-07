function fsc_options = fsc_default_options(network, model_name)

% fsc_options = fsc_default_options(network, model_name)
%
% Argument 'model_name' is optional

[nm,nr] = size(network.N);

eval(default('model_name','''model'''));

% model
fsc_options.model_name               = model_name         ; 
fsc_options.run_id                   = 'run';
fsc_options.network_CoHid            = network;

% metabolite constraints
fsc_options.fix_metabolites          = {}; % metabolites with fixed concentrations (overrides values from input model)
fsc_options.fix_metabolite_values    = [];
fsc_options.conc_min_default         = 0.001              ; % mM
fsc_options.conc_max_default         = 10                 ; % mM
fsc_options.conc_min                 = [];% 0.001 * ones(nm,1) ; % mM
fsc_options.conc_max                 = [];% 10 * ones(nm,1);   ; % mM
fsc_options.conc_fix                 = [];
fsc_options.met_fix                  = [];
fsc_options.replace_cofactors        = {};

% given data
fsc_options.c_data                   = [];
fsc_options.u_data                   = [];
fsc_options.kinetic_data             = [];

% kinetic data
fsc_options.reaction_column_names    = []; % column names (in data file) for loading of kinetic data
fsc_options.compound_column_names    = [];
fsc_options.KM_lower                 = []; % mM
fsc_options.Keq_upper                = [];
fsc_options.flag_given_kinetics      = 0;
fsc_options.kcat_usage               = 'use';
fsc_options.kcat_prior_median        = 350;   % similar to median in glycolysis+tca
fsc_options.kcat_prior_log10_std     = 0.1;   % reduce spread of kcat values
fsc_options.kcat_lower               = 50;    % 1/s
fsc_options.kcat_upper               = 2000;  % 1/s
fsc_options.GFE_fixed                = 1;     % flag
fsc_options.insert_Keq_from_data     = 0;     % flag

% specific enzyme costs
fsc_options.ind_scored_enzymes       = 1:length(network.actions);
fsc_options.enzyme_cost_weights      = ones(length(fsc_options.ind_scored_enzymes),1);
fsc_options.use_cost_weights         = 'none';

% ecm
fsc_options.initial_choice           = 'mdf'; 
fsc_options.multiple_starting_points = 0;
fsc_options.fsc_scores               = {'ecf3sp'}           ;
fsc_options.lambda_regularisation    = 10^-3; 
fsc_options.quantity_info_file       = [];
fsc_options.compute_hessian          = 0;
fsc_options.compute_elasticities     = 0;
fsc_options.compute_tolerance        = 0;
fsc_options.cost_tolerance_factor    = 1.01; % one percent
fsc_options.tolerance_from_hessian   = 0;

% graphics
fsc_options.print_graphics           = 0;
fsc_options.show_graphics            = 1;
fsc_options.show_metabolites         = network.metabolites;
