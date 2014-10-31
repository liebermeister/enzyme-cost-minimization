function fsc_save_model_and_data_sbtab(filename,network,v,r,r_orig,c_data,u_data, kinetic_data, conc_min, conc_max)

% prepare data

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

[nm,nr] = size(network.N);

formulae = network_print_formulae(network);

network.kinetics      = r;
network.kinetics.type = 'cs';
network.kinetics.c    = nan * c_data(:,1);
network.kinetics.u    = nan * u_data;


% model tables ('Compound', 'Reaction', 'QuantityData')
sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0));

% flux table
flux_table = sbtab_table_construct(struct('TableName','Flux', 'TableType','Quantity','Unit','mM/s'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','Flux'},{repmat({'flux'},nr,1),network.actions,network.reaction_KEGGID,v});

% compound GFE of formation
GFE_table = sbtab_table_construct(struct('TableName','GibbsEnergyOfFormation','TableType','Quantity','Unit','kJ/mol','StandardConcentration','1mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound','G0'},{repmat({'standard Gibbs energy of formation'},nm,1),network.metabolites, network.metabolite_KEGGID, r.mu0});

% reaction GFE table
delta_mu0      = network.N' * r.mu0;
delta_mu0_orig = kinetic_data.dmu0.median;
dGFE_table = sbtab_table_construct(struct('TableName','GibbsEnergyOfReaction','TableType','Quantity','Unit','kJ/mol','StandardConcentration','1mM'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','dG0','OriginaldG0','SumFormula'},{repmat({'standard Gibbs energy of reaction'},nr,1),network.actions, network.reaction_KEGGID, delta_mu0, delta_mu0_orig, formulae});

% metabolite concentration table
concentration_table = sbtab_table_construct(struct('TableName','Concentration', 'TableType','Quantity','Unit','mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound','Concentration'},{repmat({'concentration'},nm,1),network.metabolites, network.metabolite_KEGGID, c_data(:,1)});

% enzyme concentration table
enzyme_table = sbtab_table_construct(struct('TableName','EnzymeConcentration','TableType','Quantity','Unit','mM'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','EnzymeConcentration'},{repmat({'concentration of enzyme'},nr,1),network.actions, network.reaction_KEGGID, u_data(:,1)});

% metabolite constraint table
constraint_table = sbtab_table_construct(struct('TableName','ConcentrationConstraint', 'TableType','Quantity','Unit','mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound','ConcentrationMin','ConcentrationMax'},{repmat({'concentration'},nm,1),network.metabolites, network.metabolite_KEGGID, conc_min, conc_max});

% --------------------------------

display(sprintf('Writing model files (sbtab format) with basename\n%s', filename))

% Model files
sbtab_document_save(sbtab_document,filename,0,1);

% Data tables
sbtab_table_save(flux_table, struct('filename',[ filename '_Flux.csv']));
sbtab_table_save(concentration_table, struct('filename',[ filename '_Concentration.csv']));
sbtab_table_save(enzyme_table, struct('filename',[ filename '_EnzymeConcentration.csv']));
sbtab_table_save(GFE_table, struct('filename',[ filename '_FormationGFE.csv']));
sbtab_table_save(dGFE_table, struct('filename',[ filename '_StandardReactionGFE.csv']));
sbtab_table_save(constraint_table, struct('filename',[ filename '_ConcentrationConstraint.csv']));

% Joint table with model and all data
sbtab_document = sbtab_document_add_table(sbtab_document,'ConcentrationConstraint',constraint_table);
sbtab_document = sbtab_document_add_table(sbtab_document,'Flux',flux_table);
sbtab_document = sbtab_document_add_table(sbtab_document,'GibbsEnergyOfFormation',GFE_table);
sbtab_document = sbtab_document_add_table(sbtab_document,'GibbsEnergyOfReaction',dGFE_table);
sbtab_document = sbtab_document_add_table(sbtab_document,'Concentration',concentration_table);
sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeConcentration',enzyme_table);

sbtab_document_save_to_one(sbtab_document,[filename, '_ModelData.csv']);
