function fsc_save_model_and_data_gams(filename,network,v,r,c_data,u_data,fsc_options)

display(sprintf('Writing model files (gams format) with basename\n%s', filename))

% adjust names for sbml output

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions);

fsc_options.fix_metabolites = network_adjust_names_for_sbml_export(fsc_options.fix_metabolites);

% turn network + kinetics into sbtab format

network.kinetics = r;
network.kinetics.type = 'cs';
network_sbtab = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0));

quantity_type = network_sbtab.tables.Quantity.column.column.QuantityType;
value         = network_sbtab.tables.Quantity.column.column.Value;
compound      = network_sbtab.tables.Quantity.column.column.Compound;
reaction      = network_sbtab.tables.Quantity.column.column.Reaction;


display('Writing GAMS files');

% write table efms.csv

mytable([ [{'efms'}; network.actions], [{1}; num2cell(v)] ]','comma',[filename 'efms.csv']);

% write table metfixed.csv

mytable([column(fsc_options.fix_metabolites), num2cell(column(fsc_options.fix_metabolite_values))],'comma',[filename 'metfixed.csv']);

% write table metmoiety.csv

% write table moiety.csv

% write table newmet.csv

% write table activators.csv

ind = find(strcmp('activation constant', quantity_type));
mytable([compound(ind), reaction(ind), num2cell(value(ind))],'comma',[filename 'activators.csv']);

% write table inhibitors.csv

ind = find(strcmp('inhibition constant', quantity_type));
mytable([compound(ind), reaction(ind), num2cell(value(ind))],'comma',[filename 'inhibitors.csv']);

% write table kcats.csv

ind = find(strcmp('substrate catalytic rate constant', quantity_type));
mytable([reaction(ind), num2cell(value(ind))],'comma',[filename 'kcats.csv']);

% write table keqs.csv

ind = find(strcmp('equilibrium constant', quantity_type));
mytable([reaction(ind), num2cell(value(ind))],'comma',[filename 'keqs.csv']);

% write table kmstoich.csv

ind = find(strcmp('Michaelis constant', quantity_type));

stoich = [];
for it = 1:length(ind),
  i_met = find(strcmp(compound{ind(it)},network.metabolites));
  i_rea = find(strcmp(reaction{ind(it)},network.actions));
  stoich(it,1) = network.N(i_met,i_rea);
end

mytable([compound(ind), reaction(ind), num2cell([stoich, value(ind)])],'comma',[filename 'kmstoich.csv']);
