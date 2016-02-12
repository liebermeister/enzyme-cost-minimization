% --------------------------------------------------------
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


% --------------------------------------------------------
% Set file names


data_DIR   = [ecm_BASEDIR 'demo/data/'];

result_DIR = [ecm_BASEDIR 'demo/results/'];

infile_pb  = [data_DIR 'pb_ecoli_aerobic_minimal.tsv'];

infile_ecm = [data_DIR 'ecm_ecoli_aerobic.tsv'];
  

% --------------------------------------------------------
% Run Parameter Balancing

[report, errors] = ecm_simple(infile_pb, result_DIR, struct('actions','parameter_balancing'));


% --------------------------------------------------------
% Run Enzyme Cost Minimization

[report,errors] = ecm_simple(infile_ecm, result_DIR, struct('actions','ecm_standard'));
