% --------------------------------------------------------------------------------------
% Before running this script: run fsc_workflow_prepare_all
%
% What the script does:
% o convert the pathway models (and kinetic data from SBtab) to matlab format
% o compute GFE of formation
% o save files -> ../models/elad/[PATHWAY_COLLECTION]/matlab/

% Variables to be set: 
%  data_id            = [STRING] (identifier to be used in filenames)
%  pathway_collection = {'obd_pathways'|'kegg_modules'} 
%  organism           = {'eco'|'bsu'|'sce'}

if 0,
  data_id  = 'd1_eco'; pathway_collection = 'obd_pathways'; organism = 'eco'; fsc_workflow_prepare_single;
  data_id  = 'd2_bsu'; pathway_collection = 'obd_pathways'; organism = 'bsu'; fsc_workflow_prepare_single;
  data_id  = 'd3_sce'; pathway_collection = 'obd_pathways'; organism = 'sce'; fsc_workflow_prepare_single;
  data_id  = 'd0_eco'; pathway_collection = 'kegg_modules'; organism = 'eco'; fsc_workflow_prepare_single;
end

switch pathway_collection,
  
  case 'obd_pathways', 
    model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';

    switch organism,
      case {'eco','bsu'},
        model_names = {'EMP-GLYCOLYSIS-PTS-oxPPP-TCA'}; %'EMP-GLYCOLYSIS-PTS', 'ED-GLYCOLYSIS', 'TCA', 'EMP-GLYCOLYSIS-PTS-TCA',  'EMP-GLYCOLYSIS-PTS-TCA-OVERFLOW'};
      case 'sce',
        model_names = {'EMP-GLYCOLYSIS', 'ED-GLYCOLYSIS', 'TCA', 'EMP-GLYCOLYSIS-TCA', 'EMP-GLYCOLYSIS-PTS-TCA', 'EMP-GLYCOLYSIS-PTS-oxPPP-TCA', 'EMP-GLYCOLYSIS-TCA-OVERFLOW'};
    end
        
    %%, 'ED-GLYCOLYSIS-HIGHWAY', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'PEP-GLX-MQO', 'TCA-CHANNEL'};

    %% GEHT NICHT (singulaere matrix in enz-optimierung)
    %% 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'PEP-GLX', 'rPPP','rTCA', 'TCA-P_fluorescens' 

    %% GEHT NICHT (too restrictive constraints?)
    %%'ED-SEMI-PHOSPHORYLATIVE', 'METHYLGLYOXAL', 'oxPPP', 'PPP-PHOSPHOKETOLASE', 'rPPP-Heterotrophs',  'TCA-MQO', 'ED-NON-PHOSPHORYLATIVE', 'TCA-OA100nm', 'TCA-OA10nm', '3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'EMP-GLYCOLYSIS-HIGHWAY', 'MOP'
    
  case 'kegg_modules', 
    model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/kegg_modules/';
    model_names = {'M00001','M00002','M00003','M00004','M00005','M00006','M00007','M00008','M00009','M00010','M00011','M00012','M00013','M00014','M00015','M00016','M00017','M00018','M00019','M00020','M00021','M00022','M00023','M00024','M00025','M00026','M00027','M00028','M00029','M00030','M00031','M00032','M00033','M00034','M00035','M00036','M00037','M00038','M00039','M00040','M00042','M00043','M00044','M00045','M00046','M00047','M00048','M00050','M00051','M00052','M00053','M00055','M00060','M00061','M00063','M00064','M00066','M00067','M00081','M00082','M00083','M00086','M00087','M00088','M00089','M00090','M00091','M00092','M00093','M00094','M00095','M00096','M00097','M00098','M00099','M00100','M00101','M00102','M00103','M00104','M00106','M00107','M00108','M00109','M00110','M00112','M00113','M00114','M00115','M00116','M00117','M00118','M00119','M00120','M00121','M00122','M00123','M00124','M00125','M00126','M00127','M00128','M00129','M00130','M00131','M00132','M00133','M00134','M00135','M00136','M00137','M00138','M00139','M00140','M00141','M00142','M00144','M00145','M00148','M00149','M00150','M00151','M00152','M00165','M00166','M00167','M00168','M00169','M00170','M00171','M00172','M00173','M00174','M00175','M00176','M00307','M00308','M00309','M00310','M00311','M00312','M00313','M00338','M00344','M00346','M00347','M00356','M00357','M00358','M00364','M00365','M00366','M00367','M00368','M00369','M00370','M00371','M00372','M00373','M00374','M00375','M00376','M00377','M00378','M00415','M00416','M00418','M00419','M00422','M00432','M00433','M00525','M00526','M00527','M00528','M00529','M00530','M00531','M00532','M00533','M00534','M00535','M00537','M00538','M00539','M00540','M00541','M00543','M00544','M00545','M00547','M00548','M00551','M00555','M00561','M00562','M00563','M00567','M00568','M00569','M00570','M00572','M00573','M00577','M00579'};

end


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

  network.metabolite_names         = fsc_kegg_compound_id_to_name(network.metabolite_KEGGID, kegg_conversion_file);
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

