function ecm_save_result_sbtab(filename,network,c,u,A_forward,options)

eval(default('options','struct'));

if strcmp(filename(end-3:end),'.csv'), 
  filename(1:end-4);
end
  
options_default.flag_one_file = 1;
options_default.r             = struct;
options_default.document_name = 'Model';
options = join_struct(options_default, options);

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

[nm,nr] = size(network.N);

if ~isfield(network,'metabolite_KEGGID'),
  network.metabolite_KEGGID = network.Compound_Identifiers_kegg_compound;
end

% -------------------------------------------------------
% model with optimized metabolite and enzyme levels

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

formulae = network_print_formulae(network);

% try to insert results from enzyme prediction with common modular rate law
if isfield(c,'ecf4cmr'),
  [nm,nr] = size(network.N);
  network.kinetics      = options.r;
  network.kinetics.type = 'cs';
  network.kinetics.c    = c.ecf4cmr;
  network.kinetics.u    = u.ecf4cmr;
  sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'write_enzyme_concentrations',1,'c',c.ecf4cmr,'modular_rate_law_parameter_id', 1,'document_name', options.document_name));
else
  sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'document_name', options.document_name));
end

% model tables ('Compound', 'Reaction', 'QuantityData')
display(sprintf('Writing model file %s',  [filename '_CMR_final.csv']));
sbtab_document_save_to_one(sbtab_document,[filename '_CMR_final.csv'],0);


% -----------------------------------------------
% predicted concentrations

fn = fieldnames(c); 

c_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableType','Quantity','TableName','Predicted concentations','Document','ECM predictions','Unit','mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound'},{repmat({'concentration'},nm,1),network.metabolites, network.Compound_Identifiers_kegg_compound});

for it = 1:length(fn),
  c_table = sbtab_table_add_column(c_table,fn{it}, c.(fn{it})(:,1), 0);
end


% -----------------------------------------------
% predicted enzyme levels

fn = fieldnames(u); 

u_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableType','Quantity','TableName','Predicted enzyme levels','Document','ECM predictions','Unit','mM'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction'},{repmat({'concentration of enzyme'},nr,1),network.actions, network.Reaction_Identifiers_kegg_reaction});

for it = 1:length(fn),
  u_table = sbtab_table_add_column(u_table,fn{it}, u.(fn{it})(:,1),0);
end


% -----------------------------------------------
% predicted reaction driving forces

fn = fieldnames(A_forward); 

A_forward_table = sbtab_table_construct(struct('DocumentName', options.document_name, 'TableType','Quantity','TableName','Predicted driving forces','Document','ECM predictions','Unit','kJ/mol'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction'},{repmat({'Gibbs energy of reaction'},nr,1),network.actions, network.Reaction_Identifiers_kegg_reaction});

for it = 1:length(fn),
  A_forward_table = sbtab_table_add_column(A_forward_table,fn{it}, A_forward.(fn{it})(:,1),0);
end


% -----------------------------------------------
% save tables

switch options_default.flag_one_file, 
  
  case 0,
    display(sprintf('Writing model files (SBtab format) with basename\n%s', filename))
    sbtab_table_save(c_table, struct('filename',[ filename '_PredictedConcentrations.csv']));
    sbtab_table_save(u_table, struct('filename',[ filename '_PredictedEnzymeLevels.csv']));
    sbtab_table_save(A_forward_table, struct('filename',[ filename '_PredictedForces.csv']));
  case 1,
    display(sprintf('Writing model file (SBtab format, one file)\nBASENAME: %s', filename))
    sbtab_document = sbtab_document_construct;
    sbtab_document = sbtab_document_add_table(sbtab_document,'Concentration',c_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeConcentration',u_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'ReactionAffinity',A_forward_table);
    sbtab_document_save_to_one(sbtab_document, [ filename '_Predictions.csv']);
end
