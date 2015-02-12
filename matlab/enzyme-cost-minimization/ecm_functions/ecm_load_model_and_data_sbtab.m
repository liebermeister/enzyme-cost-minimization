function [my_network,my_v,my_c_data,my_u_data, my_conc_min, my_conc_max, my_positions,warnings] = ecm_load_model_and_data_sbtab(problem_filename_sbtab, outdir)

% [my_network,my_v,my_c_data,my_u_data, my_conc_min, my_conc_max, my_positions, warnings] = ecm_load_model_and_data_sbtab(problem_filename_sbtab)

my_sbtab = sbtab_document_load_from_one(problem_filename_sbtab);

if exist('outdir','var'),
  options.my_matlab_tmp = outdir;
end

warnings = '';

my_network   = sbtab_to_network(my_sbtab,options);

[nm,nr] = size(my_network.N);

my_v         = ones(nr,1);
my_c_data    = ones(nm,1);
my_u_data    = ones(nr,1);
my_conc_min  = ones(nm,1);
my_conc_max  = ones(nm,1);
my_positions = [];

try
  my_v         = sbtab_table_get_column(my_sbtab.tables.Flux,'Flux',1);
catch err
  warnings = 'Flux table missing';
end

try
  my_c_data    = sbtab_table_get_column(my_sbtab.tables.Concentration,'Concentration',1);
catch err
  warnings = [warnings, '; Concentration table missing'];
end

try
  my_u_data    = sbtab_table_get_column(my_sbtab.tables.EnzymeConcentration,'EnzymeConcentration',1);
catch err
  warnings = [warnings, '; Enzyme concentration table missing'];
end

try
  my_conc_min  = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'Concentration:Min',1);
  my_conc_max  = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'Concentration:Max',1);
  catch err
    warnings = [warnings, '; Concentration constraint table missing'];
end

try
  my_positions = my_sbtab.tables.Position;
catch err
  warnings = [warnings, '; Position table missing'];
end
