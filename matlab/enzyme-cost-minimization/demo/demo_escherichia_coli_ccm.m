% Demo script for Parameter Balancing and Enzyme Cost Minimization
%
% Prepared model files for E coli central metabolism.are used.
% 
% The MATLAB function 'ecm_simple' is a wrapper function 
% that reads the data, sets some default parameters, and calls 
% dedicated functions for Parameter Balancing or Enzyme Cost Minimization.
%
% The input files (in SBtab format) are given in the subdirectory 'data'
% The results (in SBtab format) are written to the subdirectory 'results'


% --------------------------------------------------------
% Set file names


data_DIR   = [ecm_BASEDIR 'demo' filesep 'data' filesep];

result_DIR = [ecm_BASEDIR 'demo' filesep 'results' filesep];


% --------------------------------------------------------
% Run Parameter Balancing

infile_pb  = [data_DIR 'ecoli_aerobic_parameter_balancing_input.tsv'];

[report, errors] = ecm_simple(infile_pb, result_DIR, struct('actions','parameter_balancing'));


% --------------------------------------------------------
% Run Enzyme Cost Minimization

infile_ecm = [data_DIR 'ecoli_aerobic_ecm_input.tsv'];

% to generate graphics, set 'make_report' to 1 instead:

[report,errors] = ecm_simple(infile_ecm, result_DIR, struct('actions','ecm_standard','make_report',1));


