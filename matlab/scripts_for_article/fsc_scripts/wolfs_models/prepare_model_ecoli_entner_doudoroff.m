% Prepare model files for the model  ecoli_entner_doudoroff_pts
%
% Save, for each, three files with: network; kinetic data; metabolic data

model_name = 'ecoli_entner_doudoroff_pts';

% MAKE SURE THAT THE DYNAMIC DATA FILE FOR THE MODEL IS UP TO DATE
% (run script_fsc3.m for the model)
% THEN RUN THE COMMANDS FOR MAKING THE KINETIC DATA (BELOW IN THIS SCRIPT, COMMENTED OUT)

make_data_files = 1;

basedir = '/home/wolfram/projekte/flux_specific_cost/matlab_flux_specific_cost/results/wolfs_models/';

    
filenames  = sbtab_model_filenames(model_name);
load(filenames.network_file);
load(filenames.flux_file);
%% scale flux to realistic value (see ~/matlab/dynamic_data/biological_numbers/glycolytic_flux)
glyc_flux_mmol_gwd_h  = 5; %1-7 mmol/(gdw * h) abb. 4 in Nanchen A, Schicker A, Sauer U (2006) Nonlinear dependency of intracellular fluxes on growth rate in miniaturized continuous cultures of Escherichia coli. Appl Environ Microbiol 72: 1164-1172 
ecoli_cell_dry_weight = 0.28 * 10^-12; % in g, from bionumbers
ecoli_cell_volume     = 1.1*10^-18; % (in m^3), from bionumbers
glyc_flux_mmol_m3_h   = glyc_flux_mmol_gwd_h  * ecoli_cell_dry_weight / ecoli_cell_volume;
glyc_flux_mmol_l_s    = 0.001 / 3600 * glyc_flux_mmol_m3_h;
v                     = v/v(1) * glyc_flux_mmol_l_s;
kinetic_data_file     = [basedir 'ecoli_entner_doudoroff_pts_my_kinetic_data'];
dynamic_data_file     = [basedir 'ecoli_entner_doudoroff_pts_my_dynamic_data'];
fsc_model_file        = [basedir '/'  model_name '_my_model'];
ind_scored_enzymes    = 1:11; % which reactions are enzyme-optimised?
conc_min              = 0.00001;
conc_max              = 100;
met_fix               = {'D_Glucose','ADP', 'ATP','Orthophosphate', 'NADH', 'NAD+', 'Acetyl_CoA','CoA','H2O'}; % 
                                                                                                               %works well
conc_fix              = [17,          0.4,  15,    8,                0.083,  2.4 ,   0.6,         1.4,  1];
%like in data, does not work
%conc_fix             = [22,           0.7,  12.2,    12.7,            0.1058,  3.24,   0.77,       1.742,  1];%];

show_metabolites  = {'D_Glucose', 'D_Glucose_6_phosphate', 'D_Glucono_1_5_lactone_6_phosphate', ...
                    '6_Phospho_D_gluconate', '2_Dehydro_3_deoxy_6_phospho_D_gluconate', ...
                    'D_Glyceraldehyde_3_phosphate', '3_Phospho_D_glyceroyl_phosphate', ...
                    '3_Phospho_D_glycerate', '2_Phospho_D_glycerate', 'Phosphoenolpyruvate', ...
                    'Pyruvate', 'Acetyl_CoA'}';

save(fsc_model_file, 'network', 'v', 'kinetic_data_file', 'dynamic_data_file', 'ind_scored_enzymes', 'conc_min', 'conc_max', 'met_fix', 'conc_fix', 'show_metabolites');

% ----------------------------------------------------------------------------------------------

if make_data_files,
  fsc_make_data_files(network, filename, kinetic_data_file, dynamic_data_file);
end
