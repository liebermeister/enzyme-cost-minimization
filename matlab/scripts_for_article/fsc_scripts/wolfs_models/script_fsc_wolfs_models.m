% --------------------------------------------------------------------------------------
% Run entire FSC analysis for several models
% --------------------------------------------------------------------------------------

% --------------------------------------------------------------------------------------
% load model and kinetic data from: 
%  script_pathway_prepare_model_arabinose_pathway.m
%  script_pathway_prepare_model_ecoli_entner_doudoroff.m
%  script_pathway_prepare_model_ecoli_glycolysis.m
 
if 0,
  clear
  % cd ~/projekte/pathway_modelling/flux_specific_cost
  % INFEASIBLE model_name = 'ecoli_glycolysis';           script_flux_specific_cost;
  model_name = 'ecoli_entner_doudoroff_pts'; script_flux_specific_cost;
  model_name = 'arabinose_pathway';          script_flux_specific_cost;
end

basedir = '/home/wolfram/projekte/flux_specific_cost/matlab_flux_specific_cost/results/wolfs_models/';

filenames.fsc_model_file    = [basedir '/' model_name '_my_model'];
filenames.kinetic_data_file = [basedir '/' model_name '_my_kinetic_data'];
filenames.dynamic_data_file = [basedir '/' model_name '_my_dynamic_data'];

fsc_scores = {'mtdf','mtdfw','mfsc2sub','fsc1','fsc2sub','fsc2','fsc3','fsc4cmr'};


% ---------------------------------------------------------
% load necessary input data (aside from fsc_scores, run_id)

load(filenames.fsc_model_file);    % network v ind_scored_enzymes met_fix 
                                   % conc_fix conc_max conc_min show_metabolites
load(filenames.kinetic_data_file); % kinetic_data
load(filenames.dynamic_data_file); % c_data u_data


% ---------------------------------------------------------
% use brenda data

run_id = 'brenda';

fsc_options.model_name         = model_name         ; 
%fsc_options.network_CoHid     = network_CoHid      ;
fsc_options.ind_scored_enzymes = ind_scored_enzymes ;
fsc_options.conc_min           = conc_min           ;
fsc_options.conc_max           = conc_max           ;
fsc_options.kinetic_data       = kinetic_data       ;
fsc_options.kcat_usage         = 'use';
%fsc_options.conc_fix           = conc_fix           ;
%fsc_options.met_fix            = met_fix            ;
fsc_options.fsc_scores         = fsc_scores         ;
fsc_options.show_metabolites   = show_metabolites   ;
fsc_options.c_data             = c_data             ;
fsc_options.u_data             = u_data             ;  
fsc_options.print_graphics     = 1                  ;  
fsc_options.run_id             = run_id             ;
fsc_options.psfile_dir = ['/home/wolfram/projekte/flux_specific_cost/ps-files/fsc_wolfs_models/' model_name '/'];

[c,u,u_tot,up,A_forward] = flux_specific_cost(network, v, fsc_options);


% ---------------------------------------------------------
% do not use brenda data

run_id                 = 'none';
fsc_options.kcat_usage = 'none';
fsc_options.run_id     = 'none';

[c,u,u_tot,up,A_forward] = flux_specific_cost(network, v, fsc_options);


% ---------------------------------------------------------
% make forward catalytic constants similar to 20

run_id                 = 'forward';
fsc_options.kcat_usage = 'forward';
fsc_options.run_id     = 'forward';

[c, u, u_tot, up, A_forward] = flux_specific_cost(network, v, fsc_options);
