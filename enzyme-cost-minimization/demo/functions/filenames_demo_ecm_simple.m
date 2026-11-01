function filenames = filenames_demo_ecm_simple()

% FILENAMES_DEMO_ECM_SIMPLE - File names and result directory for the script demo_ecm_simple
%
% filenames = filenames_demo_ecm_simple()

filenames.data_dir             = [ecm_resourcedir filesep 'models' filesep 'escherichia-coli_noor_2016' filesep];
filenames.model_data_file      = [filenames.data_dir 'escherichia-coli_noor_2016_ECM_Model.tsv'];
filenames.validation_data_file = [filenames.data_dir 'escherichia-coli_noor_2016_ECM_ValidationData.tsv'];
filenames.result_dir           = [tempdir 'demo_ecm_simple'];

[~,~] = mkdir(filenames.result_dir);
