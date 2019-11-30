function filenames = ecm_filenames()
  
filenames.demo.data_dir             = [ecm_BASEDIR filesep 'demo' filesep 'data' filesep];
filenames.demo.result_dir           = [ecm_BASEDIR filesep 'demo' filesep 'results' filesep];
filenames.demo.model_data_file      = [filenames.demo.data_dir 'ecoli_noor_2016_ECM_Model.tsv'];
filenames.demo.validation_data_file = [filenames.demo.data_dir 'ecoli_noor_2016_ECM_ValidationData.tsv'];
