%ECM_SAVE_RESULT_SBTAB - Save ECM results in SBtab format
% 
%ecm_save_result_sbtab(filename,network,c,u,A_forward,options,c_min,c_max,u_min,u_max,u_capacity,eta_energetic,eta_saturation)

function ecm_save_result_sbtab(filename,network,c,u,A_forward,options,c_min,c_max,u_min,u_max,u_capacity,eta_energetic,eta_saturation)

eval(default('options','struct','c_min', '[]', 'c_max', '[]', 'u_min', '[]', 'u_max', '[]', 'u_capacity', '[]', 'eta_energetic', '[]', 'eta_saturation', '[]'));

if strcmp(filename(end-3:end),'.tsv'), 
  filename(1:end-4);
end
  
options_default.flag_one_file = 1;
options_default.r             = struct;
options_default.document_name = 'Model';
options_default.save_tolerance_ranges = 0;
options                       = join_struct(options_default, options);

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

[nm,nr] = size(network.N);

if ~isfield(network,'metabolite_KEGGID'),
  if isfield(network,'Identifiers_kegg_compound'),
    network.metabolite_KEGGID = network.Identifiers_kegg_compound;
  end
end

if isfield(network,'metabolite_KEGGID'),
  metabolite_KEGGID = network.metabolite_KEGGID;
else
  metabolite_KEGGID = [];
end

if ~isfield(network,'reaction_KEGGID'),
  if isfield(network,'Identifiers_kegg_reaction'),
    network.reaction_KEGGID = network.Identifiers_kegg_compound;    
  end
end

if isfield(network,'reaction_KEGGID'),
  reaction_KEGGID = network.reaction_KEGGID;
else
  reaction_KEGGID = [];
end

% -------------------------------------------------------
% model with optimized metabolite and enzyme levels

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

formulae = network_print_formulae(network);

% try to insert results from enzyme prediction with common modular rate law
[nm,nr] = size(network.N);
network.kinetics      = options.r;
if isfield(c,'emc4cm'),
  network.kinetics.type = 'cs';
  network.kinetics.c    = c.emc4cm;
  network.kinetics.u    = u.emc4cm;
  sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'write_enzyme_concentrations',1,'c',c.emc4cm,'modular_rate_law_parameter_id', 1,'document_name', options.document_name));
elseif isfield(c,'emc3sp'),
  network.kinetics.type = 'ds';
  network.kinetics.c    = c.emc3sp;
  network.kinetics.u    = u.emc3sp;
  sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'write_enzyme_concentrations',1,'c',c.emc3sp,'modular_rate_law_parameter_id', 1,'document_name', options.document_name));
else
  network.kinetics.type = 'cs';
  sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'document_name', options.document_name));
end

%display(sprintf('\n'));
% model tables ('Compound', 'Reaction', 'QuantityData')
display(sprintf('Writing file %s',  [filename '_ModelState.tsv']));
sbtab_document_save_to_one(sbtab_document,[filename '_ModelState.tsv'],0);

% -----------------------------------------------
% predicted concentrations

fn = fieldnames(c); 

attributes = struct('Document', options.document_name, 'TableName','Predicted concentrations', 'TableType','Quantity', 'Unit','mM');

c_table = sbtab_table_construct(attributes, {'QuantityType','Compound'}, {repmat({'concentration'},nm,1),network.metabolites});

if length(metabolite_KEGGID),
  c_table = sbtab_table_add_column(c_table,'Compound:Identifiers:kegg.compound',metabolite_KEGGID,1);
end

for it = 1:length(fn),
  if length(c.(fn{it})),
    c_table = sbtab_table_add_column(c_table,fn{it}, c.(fn{it})(:,1), 0);
  end
end

% predicted concentration tolerance ranges
if options.save_tolerance_ranges,
  for it = 1:length(fn),
   if isfield(c_min,fn{it}),
     c_table = sbtab_table_add_column(c_table,[fn{it} ':Lower'], c_min.(fn{it})(:,1), 0);
     c_table = sbtab_table_add_column(c_table,[fn{it} ':Upper'], c_max.(fn{it})(:,1), 0);
   end
  end
end

% -----------------------------------------------
% predicted enzyme levels

fn = fieldnames(u); 

u_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableName','Predicted enzyme levels','TableType','Quantity','Document','ECM metabolic state','Unit','mM'),{'QuantityType','Reaction'},{repmat({'concentration of enzyme'},nr,1),network.actions});

if length(reaction_KEGGID),
  u_table = sbtab_table_add_column(u_table,'Reaction:Identifiers:kegg.reaction',reaction_KEGGID,1);
end

for it = 1:length(fn),
  if length(u.(fn{it})),
    u_table = sbtab_table_add_column(u_table,fn{it}, u.(fn{it})(:,1),0);
  end
end

% predicted concentration tolerance ranges
if options.save_tolerance_ranges,
  for it = 1:length(fn),
   if isfield(u_min,fn{it}),
    u_table = sbtab_table_add_column(u_table,[fn{it} ':Lower'], u_min.(fn{it})(:,1), 0);
    u_table = sbtab_table_add_column(u_table,[fn{it} ':Upper'], u_max.(fn{it})(:,1), 0);
   end
  end
end

% -----------------------------------------------
% predicted reaction driving forces

fn = fieldnames(A_forward); 

A_forward_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableName','Predicted Gibbs free energies of reaction', 'TableType','Quantity','Unit','kJ/mol'),{'QuantityType','Reaction'},{repmat({'Gibbs energy of reaction'},nr,1),network.actions});

if length(reaction_KEGGID),
  A_forward_table = sbtab_table_add_column(A_forward_table,'Reaction:Identifiers:kegg.reaction',reaction_KEGGID,1);
end

for it = 1:length(fn),
  A_forward_table = sbtab_table_add_column(A_forward_table,fn{it}, A_forward.(fn{it})(:,1),0);
end

% -----------------------------------------------
% predicted enzyme capacity (v/kcat) and efficiencies

if length(u_capacity),

  u_capacity_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableName','Predicted enzyme capacities','TableType','Quantity','Document','ECM metabolic state','Unit','kJ/mol'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','Value'},{repmat({'enzyme capacity'},nr,1),network.actions, network.reaction_KEGGID,u_capacity});
  
  fn = fieldnames(eta_energetic); 
  
  eta_energetic_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableName','Predicted energetic efficiencies','TableType','Quantity','Document','ECM metabolic state','Unit','kJ/mol'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction'},{repmat({'energetic efficiency'},nr,1),network.actions, network.reaction_KEGGID});
  
  for it = 1:length(fn),
    eta_energetic_table = sbtab_table_add_column(eta_energetic_table,fn{it}, eta_energetic.(fn{it})(:,1),0);
  end
  
  fn = fieldnames(eta_saturation); 
  
  eta_saturation_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableName','Predicted saturation efficiencies','TableType','Quantity','Document','ECM metabolic state','Unit','kJ/mol'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction'},{repmat({'saturation efficiency'},nr,1),network.actions, network.reaction_KEGGID});
  
  for it = 1:length(fn),
    eta_saturation_table = sbtab_table_add_column(eta_saturation_table,fn{it}, eta_saturation.(fn{it})(:,1),0);
  end

end

% -----------------------------------------------
% save tables

switch options_default.flag_one_file, 
  
  case 0,
    %display(sprintf('Writing model files (SBtab format) with basename\n%s', filename))
    sbtab_table_save(c_table, struct('filename',[ filename '_PredictedConcentrations.tsv']));
    sbtab_table_save(u_table, struct('filename',[ filename '_PredictedEnzymeLevels.tsv']));
    sbtab_table_save(A_forward_table, struct('filename',[ filename '_PredictedForces.tsv']));
    if length(u_capacity),
      sbtab_table_save(u_capacity_table, struct('filename',[ filename '_EnzymeCapacity.tsv']));
      sbtab_table_save(eta_energetic_table, struct('filename',[ filename '_EnergeticEfficiency.tsv']));
      sbtab_table_save(eta_saturation_table, struct('filename',[ filename '_SaturationEfficiency.tsv']));
    end
  case 1,
    %display(sprintf('Writing SBtab model file with base name: %s', filename))
    sbtab_document = sbtab_document_construct;
    sbtab_document = sbtab_document_add_table(sbtab_document,'Concentration',c_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeConcentration',u_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'ReactionAffinity',A_forward_table);
    if length(u_capacity),
      sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeCapacity',u_capacity_table);
      sbtab_document = sbtab_document_add_table(sbtab_document,'EnergeticEfficiency',eta_energetic_table);
      sbtab_document = sbtab_document_add_table(sbtab_document,'SaturationEfficiency',eta_saturation_table);
    end
    sbtab_document_save_to_one(sbtab_document, [ filename '_StateRuns.tsv']);
end
