function fsc_save_result_sbtab(filename,network,c,u,A_forward,options)

eval(default('options','struct'));

options_default.flag_one_file = 1;
options = join_struct(options_default, options);

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

[nm,nr] = size(network.N);

if ~isfield(network,'metabolite_KEGGID'),
  network.metabolite_KEGGID = network.Compound_Identifiers_kegg_compound;
end

% -----------------------------------------------
% predicted concentrations

fn = fieldnames(c); 

c_table = sbtab_table_construct(struct('TableType','Quantity','Unit','mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound'},{repmat({'concentration'},nm,1),network.metabolites, network.metabolite_KEGGID});

for it = 1:length(fn),
  c_table = sbtab_table_add_column(c_table,fn{it}, c.(fn{it})(:,1), 0);
end


% -----------------------------------------------
% predicted enzyme levels

fn = fieldnames(u); 

u_table = sbtab_table_construct(struct('TableType','Quantity','Unit','mM'),{'QuantityType','Reaction'},{repmat({'concentration of enzyme'},nr,1),network.actions});

for it = 1:length(fn),
  u_table = sbtab_table_add_column(u_table,fn{it}, u.(fn{it})(:,1),0);
end


% -----------------------------------------------
% predicted reaction driving forces

fn = fieldnames(A_forward); 

A_forward_table = sbtab_table_construct(struct('TableType','Quantity','Unit','kJ/mol'),{'QuantityType','Reaction'},{repmat({'Gibbs energy of reaction'},nr,1),network.actions});

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
    display(sprintf('Writing model file (SBtab format, one file)\n%s', filename))
    sbtab_document = sbtab_document_construct;
    sbtab_document = sbtab_document_add_table(sbtab_document,'Concentration',c_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeConcentration',u_table);
    sbtab_document = sbtab_document_add_table(sbtab_document,'ReactionAffinity',A_forward_table);
    sbtab_document_save_to_one(sbtab_document, filename);
end
