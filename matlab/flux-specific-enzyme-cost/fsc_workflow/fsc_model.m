function [c, u, u_tot, up, A_forward] = fsc_model(model_directory, data_id, run_id)

% MAKE SURE THAT THE DYNAMIC DATA FILE FOR THE MODEL IS UP TO DATE
% (run script_fsc3.m for the model)

% TEST

if 0,
  %model_directory = '/home/wolfram/projekte/flux_specific_cost/fsc/models/ecoli_glycolysis/';
  model_directory = '/home/wolfram/projekte/flux_specific_cost/fsc/models/ecoli_entner_doudoroff_pts/';
  %model_directory = '/home/wolfram/projekte/flux_specific_cost/fsc/models/arabinose_pathway/';
  data_id         = 'brenda';
  run_id          = 'run1';
  %fsc_prepare_model(model_directory);
  %fsc_prepare_data(model_directory, data_id);
  fsc_model(model_directory, data_id, run_id)
end


% -----------------------------------------------
% read prepared model and data

my_options = fsc_read_options([ model_directory '/fsc_options.txt'],fsc_file_locations);
filenames  = fsc_filenames(model_directory, model_directory, my_options.model_name, data_id,run_id);

network = [];
load(filenames.network_file);
load(filenames.flux_file);
nn = network;
load(filenames.fsc_input_file);
load(filenames.kinetic_data_file);
load(filenames.metabolic_data_file);
network = nn;


% -----------------------------------------------
% read options file 'fsc_options_run_<RUN_ID>.txt' for fsc run from model directory
% and prepare options for the FSC calculations
[ model_directory '/fsc_options_run_' run_id '.txt']
run_options = fsc_read_options([ model_directory '/fsc_options_run_' run_id '.txt']);
run_options

fsc_options.model_name         = options.model_name   ; 
fsc_options.network_CoHid      = network_CoHid;
fsc_options.ind_scored_enzymes = run_options.ind_scored_enzymes ;
fsc_options.conc_min           = run_options.conc_min           ;
fsc_options.conc_max           = run_options.conc_max           ;
fsc_options.conc_fix           = run_options.conc_fix           ;
fsc_options.met_fix            = run_options.met_fix            ;
fsc_options.fsc_scores         = run_options.fsc_scores         ;
fsc_options.show_metabolites   = run_options.show_metabolites   ;
fsc_options.c_data             = c_data;
fsc_options.u_data             = u_data;
fsc_options.print_graphics     = 1;
fsc_options.run_name           = [run_options.model_name '_' run_id];
fsc_options.kinetic_data       = kinetic_data;
fsc_options.psfile_dir         = filenames.psfile_dir;

fsc_options
% -----------------------------------------------
% Run FSC optimisation and save result (as .mat file)

[c, u, u_tot, up, A_forward, r, r_orig] = flux_specific_cost(network, v, fsc_options);

if ~length(c), return; 
end

save(filenames.fsc_result_file, 'c', 'u', 'u_tot', 'up', 'A_forward', 'r', 'r_orig', 'fsc_options');

% Save SBtab files?


% -----------------------------------------------
% Graphics

fsc_display(network,v,fsc_options,c,u,u_tot,up,A_forward,r);
