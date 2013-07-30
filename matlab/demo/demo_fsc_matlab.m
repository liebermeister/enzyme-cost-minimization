% DEMO_FSC_MATLAB - Run FSC analysis using 'flux_specific_cost.m'

fsc_dependencies;

demo_dir = [fileparts(which(mfilename))];

cd(demo_dir);

echo on;
clc
%-----------------------------------------------------------------
% Demo: Flux-specific cost calculation with predefined model
%
% Model and data files have been prepared and are read from files
% The FSC options are set within the script (using default values)
% ----------------------------------------------------------------
 
% Press key to continue
 
pause
clc
%-------------------------------------------------------------------------------
% First, we set the directory from where we load the files
% and declare the model name and ids for the set of data used 
% and the fsc run
 
model_dir  = [demo_dir '/data/ecoli_glycolysis/'];
model_name = 'ecoli_glycolysis';
data_id    = 'brenda';
run_id     = 'demo';
 
% Press key to continue
 
pause
clc
%-------------------------------------------------------------------------------
% Next, we declare the cost function(s) we are going to use
 
% Possible choices: 
% mtdf, mtdfw, fsc1, fsc2, fsc2sub, mfsc2sub, fsc3, fsc4cmr
 
% For a description, see 'help fsc_scores'
 
fsc_scores = {'fsc3'};
 
% Press key to continue
 
pause
clc
% ----------------------------------------------------------------
% We load model and kinetic data from files in the /data directory
% In the comments, we can see which variables are loaded
 
% More information on the 'network' and the 'kinetic_data'
% data structures, is given in the Metabolic Network Toolbox. 
 
filenames  = fsc_filenames(model_dir, model_dir, model_name, data_id, run_id);
 
load(filenames.network_file);         % network         Metabolic network
load(filenames.flux_file);            % v               Flux vector
load(filenames.kinetic_data_file);    % kinetic_data    Kinetic data
load(filenames.metabolic_data_file);  % c_data, u_data  Data for graphics
 
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------------
% The options for FSC are stored in a struct 'fsc_options'
% We first initialise it with default values and then insert some information
 
fsc_options = fsc_options_default(network, model_name);
 
fsc_options.c_data       = c_data;
fsc_options.u_data       = u_data;  
fsc_options.run_id       = run_id;
fsc_options.kinetic_data = kinetic_data;
 
fsc_options.conc_min     = 10^-6;        % mM
fsc_options.conc_max     = 1;            % mM
fsc_options.conc_fix     = [17, 0.4, 20, 8, 0.083, 2.4, 0.6, 1.4, 1];
fsc_options.met_fix      = {'D_Glucose','ADP', 'ATP','Orthophosphate', 'NADH', 'NAD+', 'Acetyl_CoA','CoA','H2O'}; 
 
% Press key to continue
 
pause
clc
% --------------------------------------------------------------
% Now we run the Flux-specific Cost calculation 
% with our network, flux distribution v, and options.
 
% The calculation will take a while ..
 
[c, u, u_tot, up, A_forward, r, r_orig, fsc_options] = flux_specific_cost(network, v, fsc_options);
 
% For explanations on the results, try 'help flux_specific_cost'
 
% Press key to continue
 
pause
clc
% ------------------------------------------------------------------------
% Now we plot some of the FSC results
% The command 'fsc_display' produces a number of standard plots
 
fsc_display(network,v,fsc_options,c,u,u_tot,up,A_forward,r)
 
% Press key to continue
 
pause
clc
% -------------------------------------------------------------
% That's it. All variables are still in your workspace
 
% Enjoy working with flux-specific costs!
 
% Press key to finish
pause
return
% ------------------------------------------------------------------------
% Save FSC results
 
result_dir  = model_dir;
result_file = [result_dir '/' model_name '_result_' run_id];
save(result_file, 'c', 'u', 'u_tot', 'up', 'A_forward', 'r', 'r_orig', 'fsc_options');

