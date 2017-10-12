% --------------------------------------------------------------------------------------------------
% Run ECM workflow (standard version, called by wrapper function ecm_simple) on E. coli model Noor et al 2016
%
% Use prepared model and data files 
%   data/ecoli_noor_2016_ModelData.tsv
%   data/ecoli_noor_2016_ValidationData.tsv
%
% Results are written to subdirectory "results"
%
% To see the workflow in more detail, have a look at demo_ecm_ecoli_noor_2016.m
% --------------------------------------------------------------------------------------------------

data_dir        = [ecm_BASEDIR filesep 'demo' filesep 'data' filesep];
result_dir      = [ecm_BASEDIR filesep 'demo' filesep 'results' filesep];
model_data_file = [data_dir 'ecoli_noor_2016'];

% --------------------------------------------------------
% Run Parameter Balancing with standard settings

options = struct('actions', 'parameter balancing');
[report, errors] = ecm_simple(model_data_file, result_dir, options);

% --------------------------------------------------------
% Run Enzyme Cost Minimization with standard settings

options = struct('actions', 'ecm','generate_report',1);
options.replace_cofactors = {'ATP','ADP','Orthophosphate','NADH', 'NAD+', 'NADPH','NADP+','Ubiquinone', 'Ubiquinol'};

[report, errors] = ecm_simple(model_data_file, result_dir, options);
