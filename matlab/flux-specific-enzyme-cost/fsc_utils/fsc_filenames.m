function filenames = fsc_filenames(model_directory, model_name, data_id, run_id, organism);

% filenames = fsc_filenames(model_directory, model_name, data_id, run_id, organism);

eval(default('model_name', '[]', 'data_id','[]', 'run_id','[]', 'organism','[]'));

filenames.model_dir     = model_directory;
filenames.sbtab_dir     = [model_directory '/sbtab/'];
filenames.matlab_dir    = [model_directory '/matlab/'];
filenames.data_dir      = [model_directory '/data/'];
filenames.psfile_dir    = [model_directory '/ps-files/'];

filenames.position_file = [filenames.data_dir 'network_all_Position.tsv'];
filenames.gfe_file      = [filenames.data_dir 'all_models_GFE.csv'];

filenames.resource_dir  = [fsc_root_dir '/resource'];;
filenames.kegg_conversion_file  = [filenames.resource_dir '/name-conversion/kegg_compound_names.tsv'];

if length(organism),

  filenames.protein_files  = {[filenames.data_dir 'all_models_protein_levels_' organism '.csv']};

  switch organism
    
    case 'bsu',
      filenames.protein_length_file       = [filenames.resource_dir  '/protein-properties/bsu_protein_lengths.tsv'];
      filenames.gene_name_conversion      = [filenames.resource_dir  '/name-conversion/gene_name_conversion_bsu.csv'];
      filenames.kinetic_constant_files    = {[filenames.resource_dir '/data-kinetic/kinetic_constants_brenda_bsu.tsv']};
      filenames.concentration_files       = {[filenames.resource_dir '/data-metabolite/bsu_Metabolites_BigExperiment_GM1_KEGG_first_time_point.csv']};
      filenames.enzyme_concentration_file = [filenames.resource_dir '/data-protein/Proteindata_Bsub.tsv'];
    
    case 'eco',
      filenames.protein_length_file        = [filenames.resource_dir  '/protein-properties/eco_protein_lengths.tsv'];
      filenames.gene_name_conversion       = [filenames.resource_dir  '/name-conversion/gene_name_conversion_eco.csv'];
      filenames.kinetic_constant_files     = {[filenames.resource_dir '/data-kinetic/kinetic_constants_brenda_eco.tsv']};
      filenames.concentration_files        = {[filenames.resource_dir '/data-metabolite/concentration_ecoli_glucose.tsv']; ...
                          [filenames.resource_dir '/data-metabolite/E_coli_metabolite_concentrations_Tepper_S2_v2.1_annotated.csv']; ...
                          [filenames.resource_dir '/data-metabolite/concentrations_albe_et_al_eco.csv']};
      %%filenames.enzyme_concentration_file = [filenames.resource_dir '/data-protein/E_coli_protein_mRNA.tsv'];
      filenames.enzyme_concentration_file  = [filenames.resource_dir '/data-protein/eco_abundance_table_protein_data_E_coli_Lu_1.csv'];
    
    case 'sce',
      filenames.protein_length_file       = [filenames.resource_dir  '/protein-properties/sce_protein_lengths.tsv'];
      filenames.gene_name_conversion      = [filenames.resource_dir  '/name-conversion/gene_name_conversion_sce.csv'];
      filenames.kinetic_constant_files    = {[filenames.resource_dir '/data-kinetic/kinetic_constants_brenda_sce.tsv']};
      filenames.concentration_files       = {...
                          [filenames.resource_dir '/data-metabolite/sce_oxygen_data_wiebe_aerobic_20.88.csv'];...
                          [filenames.resource_dir '/data-metabolite/sce_oxygen_data_wiebe_aerobic_2.85.csv'];...
                          [filenames.resource_dir '/data-metabolite/sce_oxygen_data_wiebe_aerobic_0.99.csv'];...
                          [filenames.resource_dir '/data-metabolite/sce_oxygen_data_wiebe_aerobic_0.5.csv'];...
                          [filenames.resource_dir '/data-metabolite/sce_oxygen_data_wiebe_anaerobic.csv'];...
                          [filenames.resource_dir '/data-metabolite/concentrations_albe_et_al_sce.csv']};
      filenames.enzyme_concentration_file = [filenames.resource_dir '/data-protein/sce_abundance_table_sce_nagaraj.csv'];
  
  end
end
  
if length(model_name),
  filenames.psfile_dir = [model_directory '/ps-files/' model_name '/'];
  filenames.model_name       = model_name;
  filenames.table_reactions  = [filenames.sbtab_dir  '/' filenames.model_name '_Reaction.tsv'];
  filenames.table_compounds  = [filenames.sbtab_dir  '/' filenames.model_name '_Compound.tsv'];
  filenames.table_positions  = [filenames.matlab_dir '/' filenames.model_name '_Position.tsv'];
  filenames.network_file     = [filenames.matlab_dir '/' filenames.model_name '_network.mat'];
  filenames.flux_file        = [filenames.matlab_dir '/' filenames.model_name '_fluxes.mat'];
  filenames.sbml_file        = [filenames.sbtab_dir  '/' filenames.model_name '_network.xml'];
  
  if length(data_id),
    filenames.kinetic_data_file   = [filenames.matlab_dir '/' filenames.model_name '_data_' data_id '_kinetic'  ];
    filenames.metabolic_data_file = [filenames.matlab_dir '/' filenames.model_name '_data_' data_id '_metabolic'];
    filenames.fsc_input_file      = [filenames.matlab_dir '/' filenames.model_name '_data_' data_id '_fsc_input'];
    
    if length(run_id),
      filenames.fsc_result_file = [filenames.matlab_dir '/' filenames.model_name '_fsc_results_' data_id '_' run_id];
    end
    
  end

end

% for graphics: 

filenames.metabolite_order_file  = [filenames.data_dir '/metabolite_order.csv'];
filenames.reaction_order_file    = [filenames.data_dir '/reaction_order.csv'];
filenames.cofactors_file         = [filenames.data_dir '/cofactors.csv']; 