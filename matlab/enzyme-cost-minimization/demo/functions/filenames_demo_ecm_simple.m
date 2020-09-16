function filenames = filenames_demo_ecm_simple()
  
filenames.data_dir             = [ecm_RESOURCEDIR filesep 'model-files' filesep 'e_coli_noor_2016' filesep];
filenames.model_data_file      = [filenames.data_dir 'e_coli_noor_2016_ECM_Model.tsv'];
filenames.validation_data_file = [filenames.data_dir 'e_coli_noor_2016_ECM_ValidationData.tsv'];
filenames.result_dir           = [tempdir 'demo_ecm_simple'];

[~,~] = mkdir(filenames.result_dir);
