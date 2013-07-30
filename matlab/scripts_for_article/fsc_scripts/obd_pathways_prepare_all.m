% ----------------------------------------------------------------
% Make joint model and use it for choosing graphics item positions 
% and Gibbs free energies of formation
%
% Compute all necessary GFE of formation and save them as SBtab file 
% ~/projekte/flux_specific_cost/matlab_flux_specific_cost/models/elad/obd_pathways/sbtab/all_models_GFE.csv
%
% (Technically, build a joint network containing all reactions, 
%  compute all GFE of formation, and save them as SBtab file)

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';

model_names = {'EMP-GLYCOLYSIS','EMP-GLYCOLYSIS-PTS', 'ED-GLYCOLYSIS', 'TCA', 'EMP-GLYCOLYSIS-TCA','EMP-GLYCOLYSIS-PTS-TCA', '3HP-4HB', '3HP', 'ACETY_COA-CHANNEL', 'DC-4HB', 'ED-GLYCOLYSIS-HIGHWAY',  'ED-SEMI-PHOSPHORYLATIVE', 'ED-NON-PHOSPHORYLATIVE','EMP-GLYCOLYSIS-HIGHWAY', 'ED-SEMI-PHOSPHORYLATIVE-HIGHWAY', 'EMP-GLYCOLYSIS-MIXED-FERMENTATION', 'METHYLGLYOXAL', 'MOP', 'oxPPP', 'PEP-GLX', 'PEP-GLX-MQO', 'PPP-PHOSPHOKETOLASE', 'rPPP', 'rPPP-Heterotrophs', 'rTCA', 'TCA-CHANNEL', 'TCA-MQO', 'TCA-OA100nm', 'TCA-OA10nm', 'TCA-P_fluorescens'};


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
network_all.metabolite_names      = kegg_compound_id_to_name(network_all.metabolite_KEGGID, kegg_conversion_file);
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

organisms = {'eco', 'bsu', 'sce'};

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

if 0,

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
  s.CompoundName = network_all.metabolites;
end
s.Value        = G0;
s.Reference    = repmat({'Component contribution method'},nm,1);
G0_sbtab = sbtab_table_construct_from_struct(struct,s,{'QuantityType','Unit','Compound MiriamID::urn.miriam.kegg:compound','Compound','Value','Reference'}');
sbtab_table_save(G0_sbtab,struct('filename',filenames.gfe_file));

else
  display('Skipping new calculation of GFE of formation');
end
