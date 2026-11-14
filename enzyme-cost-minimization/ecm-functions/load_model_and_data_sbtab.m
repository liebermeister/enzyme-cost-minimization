function [network, v, c_data, u_data, conc_min, conc_max, met_fix, conc_fix, positions, enzyme_cost_weights, network_for_plots, warnings] = load_model_and_data_sbtab(filename_model, filename_validation_data)

% LOAD_MODEL_AND_DATA_SBTAB - Load model and data from SBtab file
%
% [network, v, c_data, u_data, conc_min, conc_max, positions, warnings] = load_model_and_data_sbtab(filename_model, filename_validation_data)
%
%Load SBtab file containing (model and data) information for Enzyme Cost Minimization
%
%Data tables are matched to the model by element IDs (columns !Reaction and !Compound);
%the order of rows does not matter.
%
%For saving an SBtab file, see 'help save_model_and_data_sbtab'
%
%Arguments
% filename_model           filename for SBtab input (model and data: "..._ECM_Model.tsv"), or SBtab document object
% filename_validation_data filename for SBtab input (OPTIONAL: separate file for validation data: "..._ECM_ValidationData.tsv"), or SBtab document object
%
%Output
% network                (struct describing model, see mnt toolbox)
% v                      (nr x 1 vector of reaction rates)
% c_data                 (nm x 1 vector of measured concentrations (only for information))
% u_data                 (nr x 1 vector of measured enzyme concentrations (only for information))
% conc_min               (nm x 1 vector of minimal concentrations)
% conc_max               (nm x 1 vector of maximal concentrations)
% met_fix                (OPTIONAL: list of metabolites with fixed concentrations)
% conc_fix               (OPTIONAL: fixed concentrations corresponding to met_fix)
% positions              (OPTIONAL: positions for graphics)
% enzyme_cost_weights    (nr x 1 vector of enzyme cost weights; default [])
% network_for_plots      network struct for plots, with metabolites omitted (SBtab column !!Layout!HideInPlot) 
% warnings
  
eval(default('filename_validation_data','[]'));

my_sbtab_model      = sbtab_document_load_from_one(filename_model);

if length(filename_validation_data),
  my_sbtab_validation = sbtab_document_load_from_one(filename_validation_data);  
else
  %% no separate validation data file: look for validation data tables in the model document
  my_sbtab_validation = my_sbtab_model;
end

options = struct;

warnings = '';

network  = sbtab_to_network(my_sbtab_model,options);

[nm,nr]  = size(network.N);

v         = ones(nr,1);
c_data    = nan * ones(nm,1);
u_data    = nan * ones(nr,1);
conc_min  = nan * ones(nm,1);
conc_max  = nan * ones(nm,1);
positions = [];

enzyme_cost_weights = ones(nr,1);

reaction_id_columns = {'Reaction', 'Reaction:SBML:reaction:id'};
compound_id_columns = {'Compound', 'Compound:SBML:species:id'};

% Data tables are matched to the model elements by their IDs (columns !Reaction or !Compound);
% the order of rows does not matter. Older TableIDs are still accepted.

flux_table = get_first_table(my_sbtab_model, {'FluxData', 'Flux'});
if length(flux_table),
  v = sbtab_table_get_values_by_id(flux_table, reaction_id_columns, network.actions, 'Value');
  if any(isnan(v)), error('Flux table: fluxes missing for some reactions'); end
else
  warnings = 'Flux table missing';
end

c_table = get_first_table(my_sbtab_validation, {'MetaboliteConcentrationData', 'ConcentrationData', 'Concentration'});
if length(c_table),
  c_data = sbtab_table_get_values_by_id(c_table, compound_id_columns, network.metabolites, 'Value');
  c_data = sbtab_table_concentration_factor(c_table, 'metabolite concentration table') * c_data; % in mM
else
  warnings = [warnings, '; Concentration table missing'];
end

u_table = get_first_table(my_sbtab_validation, {'EnzymeConcentrationData', 'EnzymeData', 'EnzymeConcentration'});
if length(u_table),
  u_data = sbtab_table_get_values_by_id(u_table, reaction_id_columns, network.actions, 'Value');
  u_data = sbtab_table_concentration_factor(u_table, 'enzyme concentration table') * u_data; % in mM
else
  warnings = [warnings, '; Enzyme concentration table missing'];
end

constraint_table = get_first_table(my_sbtab_model, {'ConcentrationConstraint'});
if length(constraint_table),
  conc_min = sbtab_table_get_values_by_id(constraint_table, compound_id_columns, network.metabolites, 'Min');
  conc_max = sbtab_table_get_values_by_id(constraint_table, compound_id_columns, network.metabolites, 'Max');
  f = sbtab_table_concentration_factor(constraint_table, 'table ConcentrationConstraint');
  conc_min = f * conc_min; conc_max = f * conc_max; % in mM
else
  warnings = [warnings, '; Concentration constraint table missing'];
end

weight_table = get_first_table(my_sbtab_model, {'EnzymeCostWeight'});
if length(weight_table),
  enzyme_cost_weights = sbtab_table_get_values_by_id(weight_table, reaction_id_columns, network.actions, 'Value');
  enzyme_cost_weights(isnan(enzyme_cost_weights)) = 1;
else
  warning('Enzyme cost weight table missing in model file; I will use equal cost weights for all enzymes');
  warnings = [warnings, '; Enzyme cost weight table missing'];
end

% position table: table 'Position', or otherwise the first table of type 'Position' (e.g. TableName 'Layout')
if isfield(my_sbtab_model.tables,'Position'),
  positions = my_sbtab_model.tables.Position;
else
  table_names = fieldnames(my_sbtab_model.tables);
  for it = 1:length(table_names),
    my_attributes = my_sbtab_model.tables.(table_names{it}).attributes;
    if isfield(my_attributes,'TableType') && strcmp(my_attributes.TableType,'Position'),
      positions = my_sbtab_model.tables.(table_names{it});
      break
    end
  end
end
if isempty(positions),
  warnings = [warnings, '; Position table missing'];
end

ind      = find(conc_min == conc_max);
met_fix  = network.metabolites(ind);
conc_fix = conc_min(ind);

% network for plots

network_for_plots = network;
if length(positions),
  if sbtab_table_has_column(positions,'HideInPlot'),
    elements      = sbtab_table_get_column(positions,'Element',0);
    hide_elements = sbtab_table_get_column(positions,'HideInPlot',1);
    metabolites_to_hide = elements(find(hide_elements==1));
    network_for_plots = netgraph_simple_graph(network,metabolites_to_hide);
  end
end


% ---------------------------------------------------------------

function my_table = get_first_table(sbtab_document, table_ids)

% first table (of a list of possible TableIDs) found in an SBtab document; [] if none is found

my_table = [];
if ~isfield(sbtab_document, 'tables'), return; end
for it = 1:length(table_ids),
  if isfield(sbtab_document.tables, table_ids{it}),
    my_table = sbtab_document.tables.(table_ids{it});
    return
  end
end
