% --------------------------------------------------------------------------------------
% Run FSC analysis for a single pathway model (EMP-GLYCOLYSIS, from series 'obd_pathways')
% --------------------------------------------------------------------------------------


if 0,
  organism   = 'eco';
  model_name = 'EMP-GLYCOLYSIS-PTS'; 
  run_id = 'r1_eco_KcatBrenda';   data_id  = 'd1_eco'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r2_eco_NoKcat';       data_id  = 'd1_eco'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r3_eco_KcatForward';  data_id  = 'd1_eco'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

  organism   = 'bsu';  
  model_name = 'EMP-GLYCOLYSIS-PTS'; 
  run_id = 'r1_bsu_KcatBrenda';   data_id  = 'd2_bsu'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r2_bsu_NoKcat';       data_id  = 'd2_bsu'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r3_bsu_KcatForward';  data_id  = 'd2_bsu'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

  organism   = 'sce';  
  model_name = 'EMP-GLYCOLYSIS';
  run_id = 'r1_sce_KcatBrenda';   data_id  = 'd3_sce'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r2_sce_NoKcat';       data_id  = 'd3_sce'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
  run_id = 'r3_sce_KcatForward';  data_id  = 'd3_sce'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;
end

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';

fsc_scores = {'obd', 'obdw', 'fsc1', 'fsc2', 'fsc2sub', 'fsc3', 'fsc3prod', 'fsc4cmr'}; 

%% 'mfsc2' 
%% 'fsc4dmr' is identical to 'fsc3prod' (use just as a check)
%% 'fsc4smr' ist nicht sehr realistisch

% --------------------------------------------------------------------------------------
% load model and kinetic data (from prepared files)

filenames = fsc_filenames(model_dir, model_name, data_id, run_id, organism);

load(filenames.network_file);
load(filenames.flux_file);
load(filenames.kinetic_data_file);
load(filenames.metabolic_data_file);

v = 0.01 * v;

display(sprintf('  Rescaling the fluxes: median flux: %d mM/s',median(v)))


% ---------------------------------------------------------
% set cofactor concentrations as in data values (where these are available)

display('  Replacing cofactor concentrations');
c_data_median = nanmedian(c_data,2);
replace_cofactors =  {'ATP','ADP','Orthophosphate','NADH', 'NAD+', 'NADPH','NADP+','Ubiquinone', 'Ubiquinol'};
ll = label_names(replace_cofactors,network.metabolites);
ll = ll(find(ll));
ll = ll(find(isfinite(c_data_median(ll))));
conc_min(ll) = c_data_median(ll);
conc_max(ll) = c_data_median(ll);

% Fix glucose concentration to 12 mM

ll = label_names({'D-Glucose'},network.metabolites);
conc_min(ll) = 12;
conc_max(ll) = 12;

% Fix phosphate concentration 

ll = label_names({'Orthophosphate'},network.metabolites);
conc_min(ll) = 10;
conc_max(ll) = 10;

% Fix lactate concentration 
% ll = label_names({'(R)-Lactate'},network.metabolites);
% conc_min(ll) = 0.001;
% conc_max(ll) = 0.001;


% dabei ergibt sich noch keine realistische ATP-hydrolyse-energie!

% ---------------------------------------------------------
% prepare options 

% fsc_default options ...

clear fsc_options
fsc_options.model_name         = model_name; 
fsc_options.run_id             = run_id;
fsc_options.fsc_scores         = fsc_scores;
fsc_options.conc_min_default   = 10^-3;
fsc_options.conc_max_default   = 10^2;
fsc_options.conc_min           = conc_min;
fsc_options.conc_max           = conc_max;
fsc_options.lambda_regularisation = 10^-5;
fsc_options.kinetic_data       = kinetic_data;
fsc_options.kcat_usage         = kcat_usage;
fsc_options.kcat_prior_median  = 200; % similar to median in glycolysis+tca
fsc_options.kcat_prior_log10_std = 0.5;
fsc_options.c_data             = c_data;
fsc_options.u_data             = u_data;  
fsc_options.ind_scored_enzymes = 1:length(network.actions);
fsc_options.enzyme_cost_weights = network.enzyme_size(fsc_options.ind_scored_enzymes);
fsc_options.show_graphics      = 1;  
fsc_options.show_metabolites   = network.metabolites;
fsc_options.network_CoHid      = network_CoHid;
fsc_options.print_graphics     = 1;
fsc_options.psfile_dir         = filenames.psfile_dir;
fsc_options.quantity_info_file = [filenames.resource_dir '/data-kinetic/quantity_info.tsv'];

[c,u,u_cost,up,A_forward,r,r_orig,fsc_options] = flux_specific_enzyme_cost(network, v, fsc_options);

% check
%C = [c.fsc1, c.fsc2, c.fsc2sub, c.fsc3, c.fsc3prod, c.fsc4cmr];
%plot(C')  


save(filenames.fsc_result_file, 'model_name', 'network', 'v', 'fsc_options', 'c', 'u', 'u_cost', 'up', 'A_forward', 'r', 'r_orig', 'fsc_options');


% -------------------------------------------------------
% grafik

fsc_options.reaction_order   = load_any_table(filenames.reaction_order_file);
fsc_options.metabolite_order = load_any_table(filenames.metabolite_order_file);

fsc_display(model_name,network,v,fsc_options,c,u,u_cost,up,A_forward,r);

