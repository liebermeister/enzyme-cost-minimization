% Run standard version of the ECM workflow (wrapper function ecm_simple) on E coli model (from Noor et al 2016) 
% prepared model and data files 
%   ecoli_ccm_ProteinComposition_Haverkorn_ModelData.tsv
%   ecoli_ccm_ProteinComposition_Haverkorn_ValidationData.tsv
%
% Results are written to subdirectory "results"
%
% To see the workflow in more detail, have a look at demo_ecm_ecoli_ccm_ProteinUniform_Haverkorn.m

data_dir   = [ecm_BASEDIR filesep 'demo' filesep 'data'];
result_dir = [ecm_BASEDIR filesep 'demo' filesep 'results'];

model_data_file = [data_dir '/ecoli_ccm_ProteinComposition_Haverkorn'];

[report, errors] = ecm_simple(model_data_file,result_dir,struct('make_report',1));