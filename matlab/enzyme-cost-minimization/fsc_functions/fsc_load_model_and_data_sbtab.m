function [my_network,my_v,my_c_data,my_u_data, my_conc_min, my_conc_max] = fsc_load_model_and_data_sbtab(problem_filename_sbtab)

% [my_network,my_v,my_c_data,my_u_data, my_conc_min, my_conc_max] = fsc_load_model_and_data_sbtab(problem_filename_sbtab)

my_sbtab = sbtab_document_load_from_one(problem_filename_sbtab);

my_network  = sbtab_to_network(my_sbtab);
my_v        = sbtab_table_get_column(my_sbtab.tables.Flux,'Flux',1);
my_c_data   = sbtab_table_get_column(my_sbtab.tables.Concentration,'Concentration',1);
my_u_data   = sbtab_table_get_column(my_sbtab.tables.EnzymeConcentration,'EnzymeConcentration',1);
my_conc_min = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'MinConcentration',1);
my_conc_max = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'MaxConcentration',1);
