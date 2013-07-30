% Prepare model files for the model arabinose_pathway
%
% Save, for each, three files with: network; kinetic data; metabolic data

clear
cd ~/projekte/pathway_modelling/flux_specific_cost
model_name = 'arabinose_pathway';          flux_specific_cost_prepare_models

% MAKE SURE THAT THE DYNAMIC DATA FILE FOR THE MODEL IS UP TO DATE
% (run script_fsc3.m for the model)
% THEN RUN THE COMMANDS FOR MAKING THE KINETIC DATA (BELOW IN THIS SCRIPT, COMMENTED OUT)

make_data_files = 1;

basedir = '/home/wolfram/projekte/flux_specific_cost/matlab_flux_specific_cost/results/wolfs_models/';


filenames          = sbtab_model_filenames(model_name);
load(filenames.network_file);
load(filenames.flux_file);
    
ara_flux_mmol_l_s  = 0.35; % invented (similar to glycol. flux)
v                  = v/v(1) * ara_flux_mmol_l_s;
kinetic_data_file  = [basedir 'arabinose_pathway_my_kinetic_data'];
dynamic_data_file  = [basedir 'arabinose_pathway_my_dynamic_data'];
fsc_model_file     = [basedir '/' model_name '_my_model'];
    
ind_scored_enzymes = 1:4; % which reactions are enzyme-optimised?
conc_min           = 0.00001;
conc_max           = 100;
met_fix            = {'L_arabinose[ext]','ADP', 'ATP','Biomass'};
conc_fix           = [10, 0.4, 15, 0.001];
show_metabolites   = network.metabolites;

save(fsc_model_file, 'network', 'v', 'kinetic_data_file', 'dynamic_data_file', 'ind_scored_enzymes', 'conc_min', 'conc_max', 'met_fix', 'conc_fix', 'show_metabolites');

if make_data_files,
  fsc_make_data_files(network, filename, kinetic_data_file, dynamic_data_file);
end
