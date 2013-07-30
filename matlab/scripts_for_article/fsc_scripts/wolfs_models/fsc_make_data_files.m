function fsc_make_data_files(network, filename, kinetic_data_file, dynamic_data_file)

% fsc_make_data_files(network, filename, kinetic_data_file, dynamic_data_file)
%
% Generate data files for pathway specific cost

display('Generating data files');


% --------------------------------------------------------
%% load kinetic data 
%% brenda data originally from DATA_for_Tomer.mat  email from arren, forwarded by avi, 27.3.
  
filename = {'~/projekte/pathway_modelling/data/GFEformation_flamholz_2011_millimolar.tsv', ...
            '~/projekte/pathway_modelling/data/kinetic_constants_brenda.tsv'};

kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','Michaelis constant','activation constant',  'inhibitory constant','equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}, [], network,  filename, 0, 1);

save(kinetic_data_file,'kinetic_data');


% --------------------------------------------------------
%% load dynamic data 

%%  concentration data from table file '~/projekte/pathway_modelling/data/original_data/compound_abundance.csv' from elad
%%  protein data originally from '/home/wolfram/projekte/pathway_modelling/data/E_coli_protein_mRNA.tsv'
%%  Lu P, Vogel C, Wang R, Yao X, Marcotte EM. Absolute protein
%%  expression profiling estimates the relative contributions of
%% transcriptional and translational regulation. Nat Biotech. 2007 25(1):117-24. 

quantity_files = {'/home/wolfram/projekte/pathway_modelling/data/concentration_ecoli_Orthophosphate.tsv', ...
                  '/home/wolfram/projekte/pathway_modelling/data/concentration_ecoli_glucose.tsv'};
protein_file = {[filenames.model_dir '/' filenames.model_name '_matched_protein_mRNA.tsv']};
quantity_info_file = [];

quantity_info = data_integration_load_quantity_info({'forward enzyme mass action term','reverse enzyme mass action term','pH'}, quantity_info_file);

addpath /home/wolfram/matlab/projects/pathway_modelling/psa3
addpath /home/wolfram/matlab/projects/pathway_modelling/concentration_sampling/concentration_sampling_utils

[network, v, kinetic_data_2, quantity_info, model_quantities, basic_quantities, protein_data] = psa3_pb_prepare(filenames, [], quantity_info, quantity_files, protein_file);

c_data = kinetic_data_2.c.mean;
u_data = protein_data.u.mean;

save(dynamic_data_file,'c_data','u_data');
