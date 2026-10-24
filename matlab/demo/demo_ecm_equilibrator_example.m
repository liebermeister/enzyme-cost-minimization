% ---------------------------------------------------------------------------------------
% Enzyme cost minimization demo script
% For E coli model (Noor 2016) as implemented in the python version of ECM
%                              (available in the equilibrator python repository)  
% ---------------------------------------------------------------------------------------


% ---------------------------------------------------------------------------------------
% Set location of your ECM Model file; you can choose a different location.
% (by default, matlab will search the files on the matlab search path)

data_dir                 = [ecm_RESOURCEDIR 'model-files' filesep 'e_coli_noor_2016' filesep 'equilibrator-example' ];
filename_model           = [data_dir filesep 'e_coli_noor_2016_ecm.tsv'];
filename_validation_data = [data_dir filesep 'e_coli_noor_2016_reference.tsv'];
result_dir               = tempdir;

% ---------------------------------------------------------------------------------------
% Load model and data from the ECM Model file and translate them into matlab data structures 
% (see documentation of the Metabolic Network Toolbox for details)

display(sprintf('Reading model and data file %s', filename_model));
display(sprintf('Reading validation data file %s', filename_validation_data));

[network, v, c_data, u_data] = load_model_and_data_sbtab(filename_model, filename_validation_data);

% ---------------------------------------------------------------------------------------
% Define default options for ECM; to change the options, refer to the documentation

ecm_options            = ecm_default_options(network, 'E. coli central carbon metabolism');
ecm_options.c_data     = c_data;
ecm_options.u_data     = u_data;
ecm_options.ecm_scores = {'emc4cm'};
ecm_options            = ecm_update_options(network, ecm_options);


% ---------------------------------------------------------------------------------------
% Run ECM

[c, u, u_cost, up, A_forward, mca_info, c_min, c_max, u_min, u_max, r, u_capacity, eta_energetic, eta_saturation] = ecm_enzyme_cost_minimization(network, network.kinetics, v, ecm_options);


% ---------------------------------------------------------------------------------------
% Save results to SBtab files (and options in JSON file) (here: in your tmp directory)

document_name = 'E. coli central carbon metabolism - ECM result';

outfile_name          = [result_dir 'demo_ecoli_equilibrator'];
outfile_options_json  = [result_dir 'demo_ecoli_equilibrator_options.json'];
outfile_options_sbtab = [result_dir 'demo_ecoli_equilibrator_options.tsv' ];

options = struct('r', network.kinetics, 'method', 'emc4cm', 'document_name', document_name, 'save_tolerance_ranges', 1);

ecm_save_result_sbtab(outfile_name, network, c, u, A_forward, options, c_min, c_max, u_min, u_max, u_capacity, eta_energetic, eta_saturation, v);

% Write ECM options to JSON file
% ecm_options_save(outfile_options_json,ecm_options);

% Write ECM options to SBtab file
options_to_sbtab(ecm_options,struct('filename',outfile_options_sbtab,'TableName','Options for ECM','TableID','ConfigureECM','Method','enzyme-cost-minimization'));


% ---------------------------------------------------------------------------------------
% Display results as graphics

kinetic_data                           = [];
graphics_options.print_graphics        = 1;
graphics_options.few_graphics          = 1;
graphics_options.metabolite_order_file = [];
graphics_options.reaction_order_file   = [];
graphics_options.enzyme_colors         = sunrise_colors(length(ecm_options.ind_scored_enzymes));

ecm_display(ecm_options, graphics_options, network,v,c,u,u_cost,up,A_forward,r,kinetic_data,c_min,c_max,u_min,u_max,u_capacity,eta_energetic,eta_saturation); 
