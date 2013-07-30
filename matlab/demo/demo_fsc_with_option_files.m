% Demo script for FSC calculation, setting the options in files
%
% FSC modelling can be run by calling three functions:  
%
%   fsc_prepare_model(model_dir);          - Prepare model (from SBtab input files)
%   fsc_prepare_data(model_dir, data_id);  - Prepare kinetic data 
%   fsc_model(model_dir, run_id);          - Run FSC
%
% All files for the model are stored in a single model directory, 
% and all options are set via options files in this directory.
%
%   fsc_options.txt                  General options 
%   fsc_options_data_<DATASET>.txt   Options for a model using a specific selection of data
%   fsc_options_run_<FSC_RUN>.txt    Options for a specific FSC run
%
% For examples, see the files in our example directory

demo_dir   = [fileparts(which(mfilename))];
model_dir  = [demo_dir '/data/ecoli_glycolysis/'];
data_id    = 'brenda'; 
run_id     = 'demo'; 


fsc_prepare_model(model_dir);

% an output file is written to the model directory 


fsc_prepare_data(model_dir, data_id);

% an output file is written to the model directory 


[c, u, u_tot, up, A_forward] = fsc_model(model_dir, run_id);

% an output file is written to the model directory 
