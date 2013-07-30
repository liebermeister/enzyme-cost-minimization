% --------------------------------------------------------------------------------------
% obd_pathways_sbtab2matlab
%
% (after creating the SBtab files for Elad's models by pathway2sbtab.py:)
%
% o convert Elad's OBD pathway models and kinetic data from SBtab to matlab format
% o compute GFE of formation
% o save files -> ../models/elad/obd_pathways/matlab/

% SBtab files created by ~/projekte/flux_specific_cost/matlab_flux_specific_cost/models/elad/pathway2sbtab.py

if 0,
  data_id  = 'd1_eco'; organism = 'eco'; prepare_obd_pathways_single;
  data_id  = 'd2_bsu'; organism = 'bsu'; prepare_obd_pathways_single;
  data_id  = 'd3_sce'; organism = 'sce'; prepare_obd_pathways_single;
end

model_names = {'EMP-GLYCOLYSIS', 'EMP-GLYCOLYSIS-PTS', 'ED-GLYCOLYSIS', 'TCA', 'EMP-GLYCOLYSIS-TCA', 'EMP-GLYCOLYSIS-PTS-TCA', 'EMP-GLYCOLYSIS-TCA-OVERFLOW', 'EMP-GLYCOLYSIS-PTS-TCA-OVERFLOW'};%, 'ED-GLYCOLYSIS-HIGHWAY', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'PEP-GLX-MQO', 'TCA-CHANNEL'};

% GEHT NICHT (singulaere matrix in enz-optimierung)
% 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'PEP-GLX', 'rPPP','rTCA', 'TCA-P_fluorescens' 

% GEHT NICHT (too restrictive constraints?)
%'ED-SEMI-PHOSPHORYLATIVE', 'METHYLGLYOXAL', 'oxPPP', 'PPP-PHOSPHOKETOLASE', 'rPPP-Heterotrophs',  'TCA-MQO', 'ED-NON-PHOSPHORYLATIVE', 'TCA-OA100nm', 'TCA-OA10nm', '3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'EMP-GLYCOLYSIS-HIGHWAY', 'MOP'

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';


% ----------------------------------------------------------------
% convert models from sbtab to matlab format and save files

filenames = fsc_filenames(model_dir,[],[],[],organism);

kinetic_data_file_names    = [filenames.gfe_file; filenames.kinetic_constant_files];
protein_data_file_names    = filenames.protein_files;
metabolite_data_file_names = filenames.concentration_files;

protein_length_table = load_any_table(filenames.protein_length_file);
if strcmp(protein_length_table{1},'!!SBtab'),
  protein_length_table = protein_length_table(3:end,:);
end

kegg_conversion_file = filenames.kegg_conversion_file;

for it = 1:length(model_names),

  model_name = model_names{it};

  display(sprintf('%d: %s',it,model_name));

  filenames = fsc_filenames(model_dir, model_name, data_id,[],organism);
  
  %% convention for graphics positions: use fields
  %metabolite_KEGGID (kegg_ids) and reaction_KEGGID (sum formulae
  %with _ instead of " " 

  [network, v, conc_min, conc_max, kinetic_data] = fsc_sbtab2mnt(model_name, filenames, kinetic_data_file_names);

  network.metabolite_names         = kegg_compound_id_to_name(network.metabolite_KEGGID, kegg_conversion_file);
  network.metabolite_KEGGID        = network.metabolites;
  network.metabolites              = network.metabolite_names;
  network.graphics_par.metabolites = network.metabolite_names;
  network.graphics_par.metnames    = network.metabolite_names;
  network.reaction_formula_KEGGID          = strrep(network.formulae,' ','_');

  ind_empty = find(cellfun('length',network.metabolites)==0);
  network.metabolites(ind_empty) = network.metabolite_KEGGID(ind_empty);

  gene_name_conversion_file = filenames.gene_name_conversion;
  network.gene_names        = annotate_list(network.Gene,gene_name_conversion_file);

  for itt=1:length(network.gene_names),
    dum = strsplit('|',network.gene_names{itt});
    network.gene_names_first{itt,1} = dum{1};
  end
  ll = label_names(network.gene_names_first,protein_length_table(:,2));
  
  
  network.enzyme_size           = 350 * ones(size(network.actions));
  network.enzyme_size(find(ll)) = cell_string2num(protein_length_table(ll(find(ll)),4));
  
  network = netgraph_read_positions(network,filenames.position_file,[],1,1,network.reaction_formula_KEGGID);

  cofactors = load_any_table(filenames.cofactors_file);
  
  network_CoHid = netgraph_simple_graph(network,cofactors);
  
  gp = struct('arrowstyle','fluxes','actstyle','none','arrowvalues',v);
  figure(it);     netgraph_concentrations(network,[],v,1,gp);
  figure(100+it); netgraph_concentrations(network_CoHid,[],v,1,gp);

  %% ------------------------------------------------  
  %% protein and metabolite concentration files
  
  clear u_data c_data
  for itt = 1:length(protein_data_file_names),
    import_quantity_list = {'concentration of enzyme'};
    protein_data = data_integration_load_kinetic_data(import_quantity_list, [], network, protein_data_file_names{itt}, 0, 1,1,0);
    u_data(:,itt) = protein_data.u.median / 1000000 * 4;  
  end
  
  for itt = 1:length(metabolite_data_file_names),
    import_quantity_list = {'concentration'};
    metabolite_data = data_integration_load_kinetic_data(import_quantity_list, [], network, metabolite_data_file_names{itt}, 0, 1,1,0);
    %% assume that original data in ppm, assuming 4mM total protein 
    c_data(:,itt) = metabolite_data.c.median;
  end

  %% ------------------------------------------------  
  
  
  %% test flux stationarity
  if norm(network.N(network.external==0,:) * v),
    error('Non-stationary flux mode');
  end
  
  save(filenames.network_file,       'network', 'network_CoHid');
  save(filenames.flux_file,          'v','conc_min','conc_max');
  save(filenames.kinetic_data_file,  'kinetic_data');
  save(filenames.metabolic_data_file,'c_data','u_data');

  if ~exist(filenames.psfile_dir,'dir'),
    eval(sprintf('! mkdir %s', filenames.psfile_dir));
  end

  cd(filenames.psfile_dir);
  print(['network_all_' model_name '.eps'],'-depsc',['-f' num2str(it)]);
  print(['network_hid_' model_name '.eps'],'-depsc',['-f' num2str(100+it)]);
end

