function [network,v,c_data,u_data, conc_min, conc_max, met_fix, conc_fix,positions, enzyme_cost_weights, warnings] = ecm_load_model_and_data_sbtab(filename, tmp_dir)

% [network,v,c_data,u_data, conc_min, conc_max, positions, warnings] = ecm_load_model_and_data_sbtab(filename, tmp_dir)
%
%Load SBtab file containing (model and data) information for Enzyme Cost Minimization
%
%For saving an SBtab file, see 'help ecm_save_model_and_data_sbtab'
%
%Arguments
% filename               filename for SBtab output
% tmp_dir                a directory to which a temporary file can be written (needed for technical reasons)
%
%Output
% network                (struct describing model, see mnt toolbox)
% v                      (nr x 1 vector of reaction rates)
% r                      (struct describing model kinetics, see mnt toolbox)
% c_data                 (nm x 1 vector of measured concentrations (only for information))
% u_data                 (nr x 1 vector of measured enzyme concentrations (only for information))
% kinetic_data           (OPTIONAL: struct with kinetic data; only to give original dmu0 values)
% conc_min               (nm x 1 vector of minimal concentrations)
% conc_max               (nm x 1 vector of maximal concentrations)
% met_fix                (OPTIONAL: list of metabolites with fixed concentrations)
% conc_fix               (OPTIONAL: fixed concentrations corresponding to met_fix)
% enzyme_cost_weights    ( nr x 1 vector of enzyme cost weights; default [])
% save_single_tables     (flag for saving SBtab tables in single files; default 0)


my_sbtab = sbtab_document_load_from_one(filename);

options = struct;

if exist('outdir','var'),
  options.my_matlab_tmp = tmp_dir;
end

warnings = '';

network   = sbtab_to_network(my_sbtab,options);

[nm,nr] = size(network.N);

v         = ones(nr,1);
c_data    = nan * ones(nm,1);
u_data    = nan * ones(nr,1);
conc_min  = ones(nm,1);
conc_max  = ones(nm,1);
positions = [];
enzyme_cost_weights = ones(nr,1);

try
  v        = sbtab_table_get_column(my_sbtab.tables.Flux,'Value',1);
catch err
  warnings = 'Flux table missing';
end

try
  c_data    = sbtab_table_get_column(my_sbtab.tables.ConcentrationData,'Value',1);
catch err
  warnings = [warnings, '; Concentration table missing'];
end

try
  u_data    = sbtab_table_get_column(my_sbtab.tables.EnzymeData,'Value',1);
catch err
  warnings = [warnings, '; Enzyme concentration table missing'];
end

try
  conc_min  = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'Concentration:Min',1);
  conc_max  = sbtab_table_get_column(my_sbtab.tables.ConcentrationConstraint,'Concentration:Max',1);
  catch err
    warnings = [warnings, '; Concentration constraint table missing'];
end

try
  enzyme_cost_weight  = sbtab_table_get_column(my_sbtab.tables.EnzymeCostWeight,'Value',1);
  catch err
    warnings = [warnings, '; Enzyme cost weight table missing'];
end

try
  positions = my_sbtab.tables.Position;
catch err
  warnings = [warnings, '; Position table missing'];
end

ind = find(conc_min == conc_max);
met_fix = network.metabolites(ind);
conc_fix = conc_min(ind);
%conc_min(ind) = nan;
%conc_max(ind) = nan;