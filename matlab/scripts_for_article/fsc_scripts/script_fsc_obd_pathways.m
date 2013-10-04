% -------------------------------------------------------------------
% Run FSC for all models in model obd_pathways
% -------------------------------------------------------------------

if 0,
  model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
  data_id = 'd1_eco'; run_id = 'FSC_eco_KcatBrenda'; kcat_usage = 'use'; organism='eco';
  model_names = {'EMP-GLYCOLYSIS-PTS-oxPPP-TCA'};
  script_fsc_obd_pathways;
  %%'TCA' 'ED-GLYCOLYSIS', 'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA'

  model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
  data_id   = 'd2_bsu'; run_id = 'FSC_bsu_KcatBrenda'; kcat_usage = 'use'; organism='bsu';
  model_names = {'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA'}; % 'TCA', 
  script_fsc_obd_pathways;

  model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
  data_id   = 'd3_sce'; run_id = 'FSC_sce_KcatBrenda'; kcat_usage = 'use'; organism='sce';
  model_names = {'EMP-GLYCOLYSIS-TCA'}; %  'TCA', 'EMP-GLYCOLYSIS', 
  script_fsc_obd_pathways;
end


% GEHT AUCH: 'ED-GLYCOLYSIS-HIGHWAY', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'PEP-GLX-MQO', 'TCA-CHANNEL'

% GEHT NICHT (singulaere matrix in enz-optimierung)
% 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'PEP-GLX', 'rPPP','rTCA', 'TCA-P_fluorescens' 

% GEHT NICHT (too restrictive constraints?)
%'ED-SEMI-PHOSPHORYLATIVE', 'METHYLGLYOXAL', 'oxPPP',
%'PPP-PHOSPHOKETOLASE', 'rPPP-Heterotrophs', 'TCA-MQO', 'ED-NON-PHOSPHORYLATIVE', 
%'TCA-OA100nm', 'TCA-OA10nm', '3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'EMP-GLYCOLYSIS-HIGHWAY', 'MOP'
% geht 

fsc_scores = {'obd', 'obdw', 'fsc1', 'fsc2', 'mfsc2', 'fsc2sub','fsc3','fsc3prod','fsc4cmr'}; 

%% 'mtdfw','mfsc2sub',
%% 'fsc4dmr' is identical to 'fsc3prod' (use just as a check)
%% 'fsc4smr' nicht realistisch

for it = 1:length(model_names),

  model_name = model_names{it};

  display(sprintf('\n%d: %s',it,model_name));

  filenames = fsc_filenames(model_dir, model_name, data_id, run_id,organism);
  
  load(filenames.network_file);
  load(filenames.flux_file);
  load(filenames.kinetic_data_file);
  load(filenames.metabolic_data_file);

  %% FIX FLUXES
  display(' Using experimental flux values');
  if length(filenames.flux_data_file),  
    v_data = sbtab_table_load(filenames.flux_data_file);
    reaction_names = sbtab_table_get_column(v_data,'Reaction');
    v_values       = sbtab_table_get_column(v_data,'Flux',1);
    v  = nan * ones(size(network.actions));
    ll = label_names(network.actions,reaction_names);
    v(find(ll)) = v_values(ll(find(ll)));
    if sum(isnan(v)), v(find(isnan(v))) = nanmean(v); warning('Flux data missing'); end
  end

  v = 0.01 * v;
  display(sprintf('  Rescaling the fluxes: median flux: %d mM/s',median(v)))

  %% ---------------------------------------------------------
  %% set cofactor concentrations as in data values


  %% for which metabolites should fixed values (from model) be
  %% replaced by median data values? 

  replace_cofactors = {'ATP','ADP','Orthophosphate','NADH', 'NAD+', 'NADPH','NADP+','Ubiquinone', 'Ubiquinol'}; % ,'D-Fructose 1,6-bisphosphate', 'Citrate', 'Succinate', 'Phosphoenolpyruvate', 'D-Glucose 6-phosphate', 'D-Fructose 6-phosphate', 

  display('  Replacing cofactor concentrations');
  c_data_median = nanmean(c_data,2);

  ll = label_names(replace_cofactors,network.metabolites);
  ll = ll(find(ll));
  ll = ll(find(isfinite(c_data_median(ll))));
  conc_min(ll) = c_data_median(ll);
  conc_max(ll) = c_data_median(ll);

  %% Other fixed concentrations (mM)
  
  fix_conc = {'D-Glucose', 12; ...
              'Orthophosphate', 10; ...
              'NADH', 0.1; ...
              'NAD+', 1.2; ...
              'CoA',1; ...
              'Ubiquinol',0.5;...
              'Ubiquinone',0.1}; % 'Acetyl-CoA',0.5 };
  
  for it = 1:size(fix_conc,1),
    metab = fix_conc{it,1};
    value = fix_conc{it,2};
    ll = label_names({metab},network.metabolites);
    if ll,
      if isnan(conc_min(ll)),
        conc_min(ll) = value;
        conc_max(ll) = value;
      end
    end
  end


  %% -------------------------------------
  
  clear fsc_options
  fsc_options.model_name         = model_name; 
  fsc_options.run_id             = run_id;
  fsc_options.fsc_scores         = fsc_scores;
  fsc_options.conc_min_default   = 10^-3;
  fsc_options.conc_max_default   = 50;
  fsc_options.conc_min           = conc_min;
  fsc_options.conc_max           = conc_max;
  fsc_options.lambda_regularisation = 10^-5;
  fsc_options.kinetic_data       = kinetic_data;
  fsc_options.kcat_usage         = kcat_usage;
  fsc_options.kcat_prior_median  = 200; % similar to median in glycolysis+tca
  fsc_options.kcat_prior_log10_std = 0.1;
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

  [c,u,u_cost,up,A_forward,r,r_orig,fsc_options] = flux_specific_enzyme_cost(network,v,fsc_options);

  fsc_options.reaction_order   = load_any_table(filenames.reaction_order_file);
  fsc_options.metabolite_order = load_any_table(filenames.metabolite_order_file);

  fsc_display(model_name,network,v,fsc_options,c,u,u_cost,up,A_forward,r,kinetic_data);

  save(filenames.fsc_result_file, 'model_name', 'network', 'v', 'fsc_options', 'c', 'u', 'u_cost', 'up', 'A_forward', 'r', 'r_orig', 'fsc_options');

end

% ------------------------------------------------------------
% display:  run script_fsc_obd_pathways_display
