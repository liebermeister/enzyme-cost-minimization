% --------------------------------------------------------------------------------------------------
% Run ECM workflow on E. coli model Noor et al 2016
% (standard version of ECM, called through the wrapper function ecm_simple) 
% 
% The MATLAB function 'ecm_simple' is a wrapper function that reads the data, 
% sets default parameters, and calls functions for Parameter Balancing or Enzyme Cost Minimization.

% Use the prepared model and data files 
%   enzyme-cost-minimization/matlab/enzyme-cost-minimization/demo/data/ecoli_noor_2016_ECM_Model.tsv
%   enzyme-cost-minimization/matlab/enzyme-cost-minimization/demo/data/ecoli_noor_2016_ECM_ValidationData.tsv

% Results are written to subdirectory "enzyme-cost-minimization/matlab/enzyme-cost-minimization/demo/results"
%
% To see the workflow in more detail, have a look at demo_ecm_ecoli_noor_2016.m
% --------------------------------------------------------------------------------------------------

display(sprintf('\nThis demo script reads an existing model+data file, determines balanced model parameters, and performs ECM for a number of enzyme cost scores. Graphics are generated, and the results are saved to files.\n'))

data_dir             = [ecm_BASEDIR filesep 'demo' filesep 'data' filesep];
result_dir           = [ecm_BASEDIR filesep 'demo' filesep 'results' filesep];
model_data_file      = [data_dir 'ecoli_noor_2016_ECM_Model.tsv'];
validation_data_file = [data_dir 'ecoli_noor_2016_ECM_ValidationData.tsv'];


% --------------------------------------------------------
% Run Parameter Balancing with standard settings

options = struct('actions', 'parameter balancing','verbose',1);
[report, errors] = ecm_simple(model_data_file, validation_data_file, result_dir, options);


% --------------------------------------------------------
% Run Enzyme Cost Minimization with standard settings

options = struct('actions', 'ecm', 'generate_report',1,'verbose',1);
options.replace_cofactors = {'ATP','ADP','Orthophosphate','NADH', 'NAD+', 'NADPH','NADP+','Ubiquinone', 'Ubiquinol'};

options.fluctuations_safety_margin = 3; % safety margin (# std dev) to counter protein number fluctuations
options.cell_volume                = 1.1*10^-18;  % in m^3, default value for E coli (needed for safety margin) 

[report, errors] = ecm_simple(model_data_file, validation_data_file, result_dir, options);
