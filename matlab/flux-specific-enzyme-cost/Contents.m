% Functions for Metabolic Pathway Modelling with Flux-specific Costs
%
% MATLAB Toolboxes required: Metabolic Network Toolbox, SBML Toolbox
%
% Main functions 
%   flux_specific_cost   - FSC run for prepared model (options given in matlab) 
%   fsc_model               - Automatic FSC (1): prepare model
%   fsc_prepare_data        - Automatic FSC (2): prepare data
%   fsc_prepare_model       - Automatic FSC (3): run FSC
% 
% 
% There are two ways to run a FSC optimisation: 
% 
% 1. Create model and options within MATLAB and run 'flux_specific_cost'
% 
% 2. Prepare the model (SBtab format) and options (text files) 
%    in a model directory and run the three matlab commands 
%    'fsc_prepare_model', 'fsc_prepare_data', and 'fsc_model'.
%    All necessary information will be read from the input files.
% 
% 
% Here, option (2.) is described in more detail: 
% 
%   o Create a directory for your model ("model directory")
%   
%   o Put SBtab files describing your model into the model directory
%   
%   o Put options file 'fsc_options.txt' into the model directory
%   
%   o Put options file 'fsc_options_data_<DATASET>.txt' into the model directory
%   
%   o Put options file 'fsc_options_run_<FSC_RUN>.txt' into the directory
%   
%   o Run 'fsc_prepare_model(<MODEL_DIRECTORY>)'
%   
%   o Run 'fsc_prepare_data(<MODEL_DIRECTORY>, <DATASET>)';
%   
%   o Run 'fsc_model(<MODEL_DIRECTORY>, <FSC_RUN>)';
% 
% All output files will be written to the model directory
%    (note that files in this directory may be overwritten)
% 
% --------------------------------------------------------------------
% 
% Format of the file "fsc_options.txt"
% 
% Using the file fsc_options.txt, you can define a number of matlab variables
% that control which models will be built and how. Here is an example:
% 
% fsc_options.model_name = 'ecoli_glycolysis';
% fsc_options.sbml_file  = 'ecoli_glycolysis_network.xml';
% fsc_options.methods    = {'mtdf','mtdfw','mfsc2sub','fsc1','fsc2sub','fsc2','fsc3','fsc4cmr'};
% 
% Running the calculations:
% 
% [c, u, u_tot, up, A_forward] = fsc_model('~/projekte/flux_specific_cost/fsc/models/ecoli_glycolysis')
% Demos
%   See directory 'demo'
%
% MATLAB toolboxes required
%   Metabolic Network Toolbox (https://github.com/wolframliebermeister/mnt)
%   SBMLtoolbox    - SBML import / export  (see http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox  - SBtab format (https://github.com/wolframliebermeister/sbtab-matlab)
%   Tensor toolbox - (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
%   efmtool        - Elementary flux modes (see http://www.csb.ethz.ch/tools/efmtool)
%
% Copyright (C) 2011-2013
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>