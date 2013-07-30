% -------------------------------------------------------------------
% Run FSC for all OBD pathways
% -------------------------------------------------------------------

if 0,
  data_id   = 'd1_eco'; run_id = 'FSC_eco_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;
  data_id   = 'd2_bsu'; run_id = 'FSC_bsu_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;
  data_id   = 'd3_sce'; run_id = 'FSC_sce_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;
  %data_id   = 'd1_eco'; run_id   = 'FSC_eco_KcatForward'; kcat_usage = 'forward';
end

model_names = {'EMP-GLYCOLYSIS',  'EMP-GLYCOLYSIS-PTS', 'ED-GLYCOLYSIS', 'TCA', 'EMP-GLYCOLYSIS-TCA', 'EMP-GLYCOLYSIS-PTS-TCA'};%,  'TCA'};

% GEHT 'ED-GLYCOLYSIS-HIGHWAY', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'PEP-GLX-MQO', 'TCA-CHANNEL'

% GEHT NICHT (singulaere matrix in enz-optimierung)
% 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'PEP-GLX', 'rPPP','rTCA', 'TCA-P_fluorescens' 

% GEHT NICHT (too restrictive constraints?)
%'ED-SEMI-PHOSPHORYLATIVE', 'METHYLGLYOXAL', 'oxPPP',
%'PPP-PHOSPHOKETOLASE', 'rPPP-Heterotrophs', 'TCA-MQO', 'ED-NON-PHOSPHORYLATIVE', 
%'TCA-OA100nm', 'TCA-OA10nm', '3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'EMP-GLYCOLYSIS-HIGHWAY', 'MOP'
% geht 

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';

fsc_scores =  {'obd', 'obdw', 'fsc1', 'fsc2', 'mfsc2', 'fsc2sub','fsc3','fsc3prod','fsc4smr','fsc4cmr'}; 
%% 'mtdfw','mfsc2sub',
%% 'fsc4dmr' is identical to 'fsc3prod' (use just as a check)

for it = 1:length(model_names),

  model_name = model_names{it};

  display(sprintf('\n%d: %s',it,model_name));

  filenames = fsc_filenames(model_dir, model_name, data_id, run_id);
  
  load(filenames.network_file);
  load(filenames.flux_file);
  load(filenames.kinetic_data_file);
  load(filenames.metabolic_data_file);
  
  v = 0.01 * v;
  display(sprintf(' Rescaling the fluxes: median flux: %d mM/s',median(v)))

  %% ---------------------------------------------------------
  %% set cofactor concentrations as in data values

  display(' Replacing cofactor concentrations');
  ll = label_names({'ATP','ADP','Orthophosphate','NADH' ,'NAD+', 'NADPH','NADP+','CoA', 'Ubiquinone', 'Ubiquinol'},network.metabolites); %,
  conc_min(ll(find(ll))) = c_data(ll(find(ll)));
  conc_max(ll(find(ll))) = c_data(ll(find(ll)));

    % Fix glucose concentration to 12 mM

  ll = label_names({'D-Glucose'},network.metabolites);
  if ll,
    conc_min(ll) = 12;
    conc_max(ll) = 12;
  end

  % Fix phosphate concentration 

  ll = label_names({'Orthophosphate'},network.metabolites);
  if ll,
    conc_min(ll) = 10;
    conc_max(ll) = 10;
  end
  %% -------------------------------------
  
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
  fsc_options.kcat_standard      = 500; % similar to median in glycolysis
  fsc_options.c_data             = c_data;
  fsc_options.u_data             = u_data;
  fsc_options.ind_scored_enzymes = 1:length(network.actions);
  fsc_options.enzyme_cost_weights = network.enzyme_size(fsc_options.ind_scored_enzymes);
  fsc_options.show_graphics      = 1;  
  fsc_options.network_CoHid      = network_CoHid;
  fsc_options.show_metabolites   = network.metabolites;
  fsc_options.print_graphics     = 1;
  fsc_options.psfile_dir         = filenames.psfile_dir;
  fsc_options.quantity_info_file = [filenames.resource_dir '/data-kinetic/quantity_info.tsv'];

  [c,u,u_tot,up,A_forward,r,r_orig,fsc_options] = flux_specific_enzyme_cost(network,v,fsc_options);

  fsc_options.reaction_order   = load_any_table(filenames.reaction_order_file);
  fsc_options.metabolite_order = load_any_table(filenames.metabolite_order_file);

  fsc_display(model_name,network,v,fsc_options,c,u,u_tot,up,A_forward,r);

  save(filenames.fsc_result_file, 'model_name', 'network', 'v', 'fsc_options', 'c', 'u', 'u_tot', 'up', 'A_forward', 'r', 'r_orig', 'fsc_options');

end
