% ---------------------------------------------------------------------
% export brenda data to organism-specific data files

M = sbtab_table_load('~/projekte/flux_specific_cost/resources/kinetic_constants_brenda.tsv');

% T000..-Nummern in KEGG fuer organismen

M_bsu = sbtab_table_subselect_items(M,'OrganismIDKegg',[],'10');
M_eco = sbtab_table_subselect_items(M,'OrganismIDKegg',[],'7');
M_sce = sbtab_table_subselect_items(M,'OrganismIDKegg',[],'5');

sbtab_table_save(M_bsu,'~/projekte/flux_specific_cost/resources/bsu_kinetic_constants_brenda.tsv');
sbtab_table_save(M_eco,'~/projekte/flux_specific_cost/resources/eco_kinetic_constants_brenda.tsv');
sbtab_table_save(M_sce,'~/projekte/flux_specific_cost/resources/sce_kinetic_constants_brenda.tsv');
