
M = sbtab_table_load('~/projekte/flux_specific_cost/resources/albe_et_al_annotated_SBtab.csv');

% dann micromolar -> molar, geometrisches mittel:

conc     = 0.001 * cell_string2num(strrep(M.column.column.Concentration,'<',''));
conc_min = 0.001 * cell_string2num(M.column.column.ConcentrationMin);
conc_max = 0.001 * cell_string2num(M.column.column.ConcentrationMax);

ind = find(isnan(conc));
conc(ind) = sqrt(conc_min(ind) .* conc_max(ind));

M = sbtab_table_add_column(M,'Value',conc);
M = sbtab_table_remove_column(M,'Concentration');
M = sbtab_table_remove_column(M,'ConcentrationMin');
M = sbtab_table_remove_column(M,'ConcentrationMax');
M = sbtab_table_add_column(M,'QuantityType',repmat({'concentration'},length(conc),1));
M_sce = sbtab_table_subselect_items(M,'Organism',[],'S.cerevisiae');
M_eco = sbtab_table_subselect_items(M,'Organism',[],'E.Coli');

sbtab_table_save(M_sce,'~/projekte/flux_specific_cost/resources/sce_concentrations_albe_et_al.csv');
sbtab_table_save(M_eco,'~/projekte/flux_specific_cost/resources/eco_concentrations_albe_et_al.csv');
