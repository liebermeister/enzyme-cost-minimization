% Download the input file ecoli_ccm_ProteinUniform_Haverkorn_ModelData.tsv from 
% http://www.metabolic-economics.de/enzyme-cost-minimization/data/ecoli_ccm/ecoli_ccm_ProteinUniform_Haverkorn_ModelData.tsv

% This sets the file location of your Models and Data file; you can choose a different location.

filename = 'ecoli_noor_2016';

% This loads the model and data from the input file and translates them into
% matlab data structures (see the documentation of the Metabolic Network Toolbox for details)

[network,v,c_data,u_data, conc_min, conc_max, met_fix, conc_fix,positions, enzyme_cost_weights, warnings] = load_model_and_data_sbtab(filename);

% This defines some default options for ECM; to change the options, refer to the documentation

ecm_options = ecm_default_options(network, 'E. coli central carbon metabolism');

ecm_options.c_data = c_data;
ecm_options.u_data = u_data;

ecm_options = ecm_update_options(network, ecm_options);

% Now ECM is run

[c, u, u_cost, up, A_forward, mca_info, c_min, c_max, u_min, u_max, r, u_capacity, eta_energetic, eta_saturation] = ecm_enzyme_cost_minimization(network, network. kinetics, v, ecm_options);

% You may use this command to save all results as SBtab files (again, the file path can be changed)

document_name = 'E. coli central carbon metabolism - ECM result';
outfile_name  = 'ecoli_noor_2016_demo_ecm_results.tsv';
opt           = struct('r', network.kinetics, 'method', 'emc4cm', 'document_name', document_name, 'save_tolerance_ranges', 1);

ecm_save_result_sbtab(outfile_name, network, c, u, A_forward, opt, c_min, c_max, u_min, u_max, u_capacity, eta_energetic, eta_saturation);

% To display graphical output, use the following lines:

kinetic_data                           = [];
ecm_options.show_graphics              = 1;
graphics_options.print_graphics        = 1;
graphics_options.few_graphics          = 1;
graphics_options.metabolite_order_file = [];
graphics_options.reaction_order_file   = [];
graphics_options.enzyme_colors         = sunrise_colors(length(ecm_options.ind_scored_enzymes));
ecm_display(ecm_options, graphics_options, network,v,c,u,u_cost,up,A_forward,r,kinetic_data,c_min,c_max,u_min,u_max,u_capacity,eta_energetic,eta_saturation); 