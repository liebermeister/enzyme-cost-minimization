function fsc_options = fsc_options_default(network, model_name)

%   Argument 'model_name' is optional
%
%   fsc_options.model_name:           'model'
%   fsc_options.run_id:               'run'
%   fsc_options.network_CoHid:        network
%   fsc_options.fsc_scores:           {'fsc2'}
%   fsc_options.ind_scored_enzymes = 1:length(network.actions);
%   fsc_options.conc_min:             10^-6
%   fsc_options.conc_max:             10^2
%   fsc_options.conc_fix:             []
%   fsc_options.met_fix:              []
%   fsc_options.c_data:               []
%   fsc_options.u_data:               []
%   fsc_options.kinetic_data:         []
%   fsc_options.show_metabolites      network.metabolites;
%   fsc_options.print_graphics:       0
%   fsc_options.flag_given_kinetics:  0
%   fsc_options.enzyme_cost_weights: ones(size(ind_scored_enzymes))

eval(default('model_name','''model'''));

fsc_options.model_name         = model_name         ; 
fsc_options.run_id             = 'run'              ;
fsc_options.network_CoHid      = network;
fsc_options.fsc_scores         = {'ecf3sp'}           ;
fsc_options.ind_scored_enzymes = 1:length(network.actions);
fsc_options.conc_min_default   = 0.001              ; % mM
fsc_options.conc_max_default   = 10                 ; % mM
fsc_options.conc_min           = 0.001              ; % mM
fsc_options.conc_max           = 1                  ; % mM
fsc_options.conc_fix           = []                 ;
fsc_options.met_fix            = []                 ;
fsc_options.c_data             = []                 ;
fsc_options.u_data             = []                 ;  
fsc_options.show_metabolites   = network.metabolites;
fsc_options.print_graphics     = 0                  ;  
fsc_options.kinetic_data       = [];                ;
fsc_options.flag_given_kinetics =  0;
fsc_options.enzyme_cost_weights = ones(length(fsc_options.ind_scored_enzymes),1);
fsc_options.replace_cofactors  = {};
fsc_options.use_cost_weights  = 'none';
fsc_options.insert_original_equilibrium_constants = 0;
fsc_options.variation_for_relaxed_optimality = 0;
fsc_options.variability_u_factor             = 1.025;  