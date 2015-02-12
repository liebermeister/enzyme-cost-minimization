function ecm_save_model_and_data_sbtab(filename,network,v,r,c_data,u_data, kinetic_data, conc_min, conc_max,document_name)


eval(default('v','[]','r','[]','c_data','[]','u_data','[]','kinetic_data','[]','conc_min','[]','conc_max','[]'));

% -------------------------------------------------------
% prepare data

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

[nm,nr] = size(network.N);

formulae = network_print_formulae(network);

network.kinetics      = r;
network.kinetics.type = 'cs';
network.kinetics.c    = nan * c_data(:,1);
network.kinetics.u    = nan * u_data;

% model tables ('Compound', 'Reaction', 'QuantityData')
sbtab_document = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'document_name',document_name));
sbtab_document_save(sbtab_document,filename,0,1);

% flux table
if length(v),      
  flux_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','Flux', 'TableType','Quantity','Unit','mM/s'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','Flux'},{repmat({'flux'},nr,1),network.actions,network.reaction_KEGGID,v});
  sbtab_table_save(flux_table, struct('filename',[ filename '_Flux.csv'])); 
  sbtab_document = sbtab_document_add_table(sbtab_document,'Flux',flux_table);
end

% compound GFE of formation
if length(r),
  GFE_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','GibbsEnergyOfFormation','TableType','Quantity','Unit','kJ/mol','StandardConcentration','1mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound','Value'},{repmat({'standard Gibbs energy of formation'},nm,1),network.metabolites, network.metabolite_KEGGID, r.mu0});
  sbtab_table_save(GFE_table, struct('filename',[ filename '_FormationGFE.csv'])); 
  sbtab_document = sbtab_document_add_table(sbtab_document,'GibbsEnergyOfFormation',GFE_table);
end

% reaction GFE table
if length(r),
  delta_mu0      = network.N' * r.mu0;
  if length(kinetic_data),
    delta_mu0_orig = kinetic_data.dmu0.median;
  else
    delta_mu0_orig = nan*delta_mu0;
  end
  dGFE_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','GibbsEnergyOfReaction','TableType','Quantity','Unit','kJ/mol','StandardConcentration','1mM'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','Value','OriginalValue','SumFormula'},{repmat({'standard Gibbs energy of reaction'},nr,1),network.actions, network.reaction_KEGGID, delta_mu0, delta_mu0_orig, formulae});
  sbtab_table_save(dGFE_table, struct('filename',[ filename '_StandardReactionGFE.csv'])); 
  sbtab_document = sbtab_document_add_table(sbtab_document,'GibbsEnergyOfReaction',dGFE_table);
end

% metabolite concentration table
if length(c_data), 
  concentration_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','Concentration', 'TableType','Quantity','Unit','mM'),...
                                              {'QuantityType','Compound','Compound:Identifiers:kegg.compound','Concentration'},...
                                              {repmat({'concentration'},nm,1),network.metabolites, network.metabolite_KEGGID, c_data(:,1)});
  sbtab_table_save(concentration_table, struct('filename',[ filename '_Concentration.csv'])); 
  sbtab_document = sbtab_document_add_table(sbtab_document,'Concentration',concentration_table);
end

% enzyme concentration table
if length(u_data),
  enzyme_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','EnzymeConcentration','TableType','Quantity','Unit','mM'),{'QuantityType','Reaction','Reaction:Identifiers:kegg.reaction','EnzymeConcentration'},{repmat({'concentration of enzyme'},nr,1),network.actions, network.reaction_KEGGID, u_data(:,1)});
  sbtab_table_save(enzyme_table, struct('filename',[ filename '_EnzymeConcentration.csv']));  
  sbtab_document = sbtab_document_add_table(sbtab_document,'EnzymeConcentration',enzyme_table);
end

% metabolite constraint table
if length(conc_min),
  constraint_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','ConcentrationConstraint', 'TableType','Quantity','Unit','mM'),{'QuantityType','Compound','Compound:Identifiers:kegg.compound','Concentration:Min','Concentration:Max'},{repmat({'concentration'},nm,1),network.metabolites, network.metabolite_KEGGID, conc_min, conc_max});
  sbtab_table_save(constraint_table, struct('filename',[ filename '_ConcentrationConstraint.csv']));
  sbtab_document = sbtab_document_add_table(sbtab_document,'ConcentrationConstraint',constraint_table);
end

% position table
x = network.graphics_par.x(1,:)';
y = network.graphics_par.x(2,:)';
position_table = sbtab_table_construct(struct('DocumentName',document_name,'TableName','Position', 'TableType','Position'),{'Element','PositionX','PositionY'},{[network.metabolites; network.actions],x,y});
sbtab_table_save(position_table, struct('filename',[ filename '_Position.csv']));
sbtab_document = sbtab_document_add_table(sbtab_document,'Position',position_table);


% --------------------------------

sbtab_document_save_to_one(sbtab_document,[filename, '_ModelData.csv']);

display(sprintf('Wrote model files (sbtab format) with basename\n%s', filename))

