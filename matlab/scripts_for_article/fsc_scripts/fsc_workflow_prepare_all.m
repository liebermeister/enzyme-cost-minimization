% ----------------------------------------------------------------
% Script for preparing the ECF calculation: 
%
%   o Make joint model
%   o Use it for choosing graphics item positions 
%   o Use it for finding the Gibbs free energies of formation
%     (Compute all necessary GFE of formation and save them as SBtab file in)
%      ~/projekte/flux_specific_cost/matlab_flux_specific_cost/models/elad/[PATHWAY_SET]/sbtab/all_models_GFE.csv
%     (Technically, build a joint network containing all reactions, 
%      compute all GFE of formation, and save them as SBtab file)
%
% Variables to be set: 
%  pathway_collection = {'obd_pathways'|'kegg_modules'} 
%  organism           = {'eco'|'bsu'|'sce'}
%  flag_thermo        = 0|1   (compute thermodynamic data?)
%
% Before running this script: 
%  translate pathways from Elad's format into SBtab format
%  (see /home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/README)
%  python script pathway2sbtab.py - > files in ../models/elad/[PATHWAY_COLLECTION]/sbtab/


if 0,
  pathway_collection = 'obd_pathways'; organisms = {'eco', 'bsu', 'sce'}; flag_thermo = 0; fsc_workflow_prepare_all;
  pathway_collection = 'kegg_modules'; organisms = {'eco'};               flag_thermo = 0; fsc_workflow_prepare_all;
end


  
% ---------------------------------------------------------

switch pathway_collection,
  
  case 'obd_pathways',
    model_dir   = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
    model_names = {'3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'DIHYDROXYACETONE', 'ED-GLYCOLYSIS', 'ED-GLYCOLYSIS-HIGHWAY', 'ED-NON-PHOSPHORYLATIVE', 'ED-SEMI-PHOSPHORYLATIVE', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'EMP-GLYCOLYSIS', 'EMP-GLYCOLYSIS-HIGHWAY', 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA', 'EMP-GLYCOLYSIS-PTS-oxPPP-TCA', 'EMP-GLYCOLYSIS-PTS-TCA-OVERFLOW', 'EMP-GLYCOLYSIS-TCA', 'EMP-GLYCOLYSIS-TCA-OVERFLOW', 'METHYLGLYOXAL', 'MOP', 'oxPPP', 'PEP-GLX', 'PEP-GLX-MQO', 'PPP-PHOSPHOKETOLASE', 'REDUCTIVE_GLYCINE', 'RIBULOSE_MPP', 'rPPP', 'rPPP-Heterotrophs', 'rTCA', 'SERINE', 'SERINE-PYRUVATE', 'TCA-CHANNEL', 'TCA', 'TCA-MQO', 'TCA-OA100nm', 'TCA-OA10nm', 'TCA-P_fluorescens'}';

  case 'kegg_modules',
    model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/kegg_modules/';
    model_names = {'M00001','M00002','M00003','M00004','M00005','M00006','M00007','M00008','M00009','M00010','M00011','M00012','M00013','M00014','M00015','M00016','M00017','M00018','M00019','M00020','M00021','M00022','M00023','M00024','M00025','M00026','M00027','M00028','M00029','M00030','M00031','M00032','M00033','M00034','M00035','M00036','M00037','M00038','M00039','M00040','M00042','M00043','M00044','M00045','M00046','M00047','M00048','M00050','M00051','M00052','M00053','M00055','M00060','M00061','M00063','M00064','M00066','M00067','M00081','M00082','M00083','M00086','M00087','M00088','M00089','M00090','M00091','M00092','M00093','M00094','M00095','M00096','M00097','M00098','M00099','M00100','M00101','M00102','M00103','M00104','M00106','M00107','M00108','M00109','M00110','M00112','M00113','M00114','M00115','M00116','M00117','M00118','M00119','M00120','M00121','M00122','M00123','M00124','M00125','M00126','M00127','M00128','M00129','M00130','M00131','M00132','M00133','M00134','M00135','M00136','M00137','M00138','M00139','M00140','M00141','M00142','M00144','M00145','M00148','M00149','M00150','M00151','M00152','M00165','M00166','M00167','M00168','M00169','M00170','M00171','M00172','M00173','M00174','M00175','M00176','M00307','M00308','M00309','M00310','M00311','M00312','M00313','M00338','M00344','M00346','M00347','M00356','M00357','M00358','M00364','M00365','M00366','M00367','M00368','M00369','M00370','M00371','M00372','M00373','M00374','M00375','M00376','M00377','M00378','M00415','M00416','M00418','M00419','M00422','M00432','M00433','M00525','M00526','M00527','M00528','M00529','M00530','M00531','M00532','M00533','M00534','M00535','M00537','M00538','M00539','M00540','M00541','M00543','M00544','M00545','M00547','M00548','M00551','M00555','M00561','M00562','M00563','M00567','M00568','M00569','M00570','M00572','M00573','M00577','M00579'};

end


% ------------------------------------------------------------
% combine all networks -> network_all

filenames = fsc_filenames(model_dir);

all_reactions         = {};
all_KEGG_reaction_ids = {};

network_all = sbtab_to_network([filenames.sbtab_dir '/' model_names{1}], struct('load_quantity_table',0));

formulae        = {};
reaction_KEGGID = {};
Genes           = {};
reaction_list   = {};

for it1 = 1:length(model_names),
  display(sprintf('%d/%d',it1,length(model_names)))
  my_network      = sbtab_to_network([filenames.sbtab_dir '/' model_names{it1}], struct('load_quantity_table',0));
  formulae        = [formulae;        my_network.formulae];
  reaction_KEGGID = [reaction_KEGGID; my_network.MiriamID__urn_miriam_kegg_reaction];
  Genes           = [Genes;           my_network.Gene];
  reaction_list   = [reaction_list;   my_network.actions];
end

rr              = unique(formulae);
ind             = label_names(rr,formulae);
reaction_list   = reaction_list(ind);
Genes           = Genes(ind);
formulae        = formulae(ind);
reaction_KEGGID = reaction_KEGGID(ind);
KEGG_formula_ID = strrep(formulae,' ','_');

kegg_conversion_file = filenames.kegg_conversion_file;

network_all = network_build_from_sum_formulae_list(formulae);
network_all.actions               = reaction_list;
network_all.Gene                  = Genes;
network_all.metabolite_KEGGID     = network_all.metabolites;
network_all.metabolite_names      = fsc_kegg_compound_id_to_name(network_all.metabolite_KEGGID, kegg_conversion_file);
network_all.reaction_KEGGID       = reaction_KEGGID;
network_all.reaction_formula_KEGGID = KEGG_formula_ID;
network_all.metabolites           = network_all.metabolite_names;
network_all.graphics_par.metnames = network_all.metabolite_names;
network_all.graphics_par.actnames = network_all.actions;


% ------------------------------------------------------------
%% convention for graphics positions: use fields
%metabolite_KEGGID (kegg_ids) and reaction_formula_KEGGID (sum formulae
%with _ instead of " " 

cofactors = load_any_table(filenames.cofactors_file);

network_all_CoHid = netgraph_simple_graph(network_all,cofactors);

% figure(1000); netgraph_edit_positions(network_all,filenames.position_file,1,struct,network_all.reaction_formula_KEGGID);

% figure(1000); netgraph_edit_positions(network_all_CoHid,filenames.position_file,1,struct,network_all.reaction_formula_KEGGID);

save([filenames.data_dir '/network_all'],'network_all', 'network_all_CoHid', 'cofactors');


% ------------------------------------------------------------
% extract protein level data and save to file (separate for each organism)

for it_org = 1:length(organisms),
  
organism = organisms{it_org};

filenames = fsc_filenames(model_dir,[],[],[],organism);

gene_name_conversion_file = filenames.gene_name_conversion;

network_all.gene_names = annotate_list(network_all.Gene,gene_name_conversion_file);
network_all.gene_names(find(strcmp('-',network_all.gene_names))) = network_all.Gene(find(strcmp('-',network_all.gene_names)));
network_all.gene_names(find(strcmp('',network_all.gene_names))) = network_all.Gene(find(strcmp('',network_all.gene_names)));

% take care of multiple gene names (separated by '|')
% For each reaction, genes should be ordered according to their importance 
% (important to less important); here when making the protein data file,
% only the first  gene is used AND WILL BE THE ONLY ONE
% TO BE READ WHEN THE FILE IS PARSED to bild the kinetic_data structure.
% AN ALTERNATIVE WOULD BE TO ADD THE VALUES FOR ALL GENES BELONGING TO A REACTION ALREADY HERE

all_gene_names = {};
all_reactions  = {};
all_reactions_KEGGID = {};

for it = 1:length(network_all.gene_names),
  mm = strsplit('|',network_all.gene_names{it});
  all_gene_names = [all_gene_names; mm(1)];
  all_reactions  = [all_reactions; network_all.actions(it)];
  all_reactions_KEGGID  = [all_reactions_KEGGID; network_all.reaction_KEGGID(it)];
end

t = sbtab_table_load(filenames.enzyme_concentration_file);
t = sbtab_table_add_annotation(t,'Gene','Reaction',all_gene_names,all_reactions);
t = sbtab_table_add_annotation(t,'Gene','Reaction MiriamID::urn.miriam.kegg:reaction',all_gene_names,all_reactions_KEGGID);
t = sbtab_table_subselect_items(t,'Reaction');

sbtab_table_save(t,struct('filename',filenames.protein_files{1}));

end


% ------------------------------------------------------------
% collect GFE of formation

if flag_thermo,

  filenames = fsc_filenames(model_dir);  load([filenames.data_dir '/network_all']);

  G0 = network_gfe_of_formation(network_all);

  %% [G0,GO_sbtab] = network_gfe_of_formation(network_all);
  %% FUNKTIONIERT NUR; WENN ES DIREKT AUSGEFÜHRT WIRD: 
  %% (NICHT INNERHALB network_gfe_of_formation)
  
  nm = length(network_all.metabolites);
  clear s
  s.QuantityType = repmat({'standard chemical potential'},nm,1);
  s.Unit         = repmat({'kJ/mol'},nm,1);
  s.Compound_MiriamID__urn_miriam_kegg_compound = network_all.metabolite_KEGGID;
  %s.SBMLSpeciesID = network_all.metabolite_KEGGID;
  if isfield(network_all,'metabolite_names'),
    s.CompoundName = network_all.metabolite_names;
  else,
    s.CompoundName = network_all.metabolites; % hier: KEGG2ID
  end
  s.Value        = G0;
  s.Reference = repmat({'Component contribution method'},nm,1);
  G0_sbtab    = sbtab_table_construct_from_struct(struct,s,{'QuantityType','Unit','Compound MiriamID::urn.miriam.kegg:compound','Compound','Value','Reference'}');
  sbtab_table_save(G0_sbtab,struct('filename',filenames.gfe_file));
  
else
  display('No new calculation of GFE of formation');
end

