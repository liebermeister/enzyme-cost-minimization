function ecm_save_model_and_data_gams(filename,network,v,r,c_data,u_data,enzyme_cost_weights,ecm_options)

% ecm_save_model_and_data_gams(filename,network,v,r,c_data,u_data,enzyme_cost_weights,ecm_options)
% 
% Convert data for Enzyme Cost Minimization (model and data) from SBtab format to GAMS input format
% 
% For generating the input file (SBtab format), see 'help ecm_save_model_and_data_sbtab'

display(sprintf('Writing input files for the NEOS optimization server to directory\n%s', filename))

% adjust names for sbml output

[network.metabolites,network.actions] = network_adjust_names_for_sbml_export(network.metabolites,network.actions,1);

% ecm_options.fix_metabolites = network_adjust_names_for_sbml_export(ecm_options.fix_metabolites,[],1);

% turn network + kinetics into sbtab format

network.kinetics = r;
network.kinetics.type = 'cs';
network_sbtab = network_to_sbtab(network, struct('use_sbml_ids',0,'verbose',0,'write_concentrations',0,'write_enzyme_concentrations',0));

quantity_type = network_sbtab.tables.Quantity.column.column.QuantityType;
value         = network_sbtab.tables.Quantity.column.column.Value;
compound      = network_sbtab.tables.Quantity.column.column.Compound;
reaction      = network_sbtab.tables.Quantity.column.column.Reaction;

[compound,reaction] = network_adjust_names_for_sbml_export(compound,reaction,1);
%display('Writing GAMS files');

% write table efms.csv (DEPRECATED)
% mytable([ [{'efms'}; network.actions], [{1}; num2cell(v)] ]','comma',[filename 'efms.csv']);

% write table metfixed.csv

mytable([column(ecm_options.met_fix), num2cell(column(ecm_options.conc_fix))],'comma',[filename 'metfixed.csv']);

% write table boundse.csv

[nm,nr] = size(network.N);

mytable([[network.actions, repmat({'lo','0'},nr,1)]; ...
         [network.actions, repmat({'up','100'},nr,1)]],'comma',[filename 'boundse.csv']);

% write table boundsx.csv

mytable([[network.metabolites, repmat({'lo'},nm,1), num2cell(ecm_options.conc_min)]; ...
         [network.metabolites, repmat({'up'},nm,1), num2cell(ecm_options.conc_max)]],'comma',[filename 'boundsx.csv']);

% write table moietymet.csv

%mytable({[]},'comma',[filename 'moietymet.csv']);

% write table moietyval.csv

%mytable({[]},'comma',[filename 'moietyval.csv']);

% write table alpha.csv

%mytable({[]},'comma',[filename 'alpha.csv']);

% write table a.csv

%mytable({[]},'comma',[filename 'a.csv']);

% write table newmet.csv (DEPRECATED)
% mytable({[]},'comma',[filename 'newmet.csv']);

% write table activators.csv

%ind = find(strcmp('activation constant', quantity_type));
%mytable([compound(ind), reaction(ind), num2cell(value(ind))],'comma',[filename 'activators.csv']);

% write table inhibitors.csv

%ind = find(strcmp('inhibition constant', quantity_type));
%mytable([compound(ind), reaction(ind), num2cell(value(ind))],'comma',[filename 'inhibitors.csv']);

% write table rates.csv

mytable([[{'0'}, network.actions']; [{'1'}, num2cell(v)']],'comma',[filename 'rates.csv']);

% write table weights.csv

mytable([network.actions, num2cell(enzyme_cost_weights)],'comma',[filename 'weights.csv']);

% write table kcats.csv

ind = find(strcmp('substrate catalytic rate constant', quantity_type));
mytable([reaction(ind), num2cell(value(ind))],'comma',[filename 'kcats.csv']);

% write table keqs.csv

ind = find(strcmp('equilibrium constant', quantity_type));
mytable([reaction(ind), num2cell(value(ind))],'comma',[filename 'keqs.csv']);

% write tables stoich.csv and kms.csv

ind = find(strcmp('Michaelis constant', quantity_type));

stoich = [];
for it = 1:length(ind),
  i_met = find(strcmp(compound{ind(it)},network.metabolites));
  i_rea = find(strcmp(reaction{ind(it)},network.actions));
  stoich(it,1) = network.N(i_met,i_rea);
end

mytable([reaction(ind), compound(ind), num2cell([stoich])],'comma',[filename 'stoich.csv']);

mytable([reaction(ind), compound(ind), num2cell([value(ind)])],'comma',[filename 'kms.csv']);

%eval(sprintf(' cd %s; zip MODEL.zip *tsv;', filename));

display(sprintf('To run an optimisation with these files, please use the NEOS server at\nhttps://proto.neos-server.org/neos/solvers/application:MER/csv.html'))


