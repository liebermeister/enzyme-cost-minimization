function [network, v, conc_min, conc_max, kinetic_data] = fsc_sbtab2mnt(model_name, filenames,kinetic_data_file_names,options,organism_long)

% [network, v, conc_fix, kinetic_data] = fsc_sbtab2mnt(model_name, model_dir, matlab_dir,kinetic_data_file_names)
%
% Read pathways (SBtab format), collect dG0' values (from network_thermo)
% and convert them into MNT network format

% I use mM and convert everything from Elad's convention within this script

% if organism_long is given, the kinetic data are filtered for the organism

eval(default('options','struct'));
options_default = struct();
options = join_struct(options_default,options);

network            = sbtab_to_network([filenames.sbtab_dir '/' model_name], struct('load_quantity_table',0));
v                  = cell_string2num(network.Flux);
network.external   = [network.N * v ~= 0];
network.metabolite_names = fsc_kegg_compound_id_to_name(network.metabolite_KEGGID,filenames.kegg_conversion_file);


% -----------------------------------------------------------
% omit protons from network model

ind_proton = find(strcmp('C00080',network.metabolite_KEGGID));

if ind_proton,
  network = network_remove(network,ind_proton);
end


% -----------------------------------------------------------
% extract kegg reaction IDs

for it=1:length(network.actions),
   dum = strsplit('_',network.actions{it});
   network.reaction_KEGGID{it,1} = dum{2};
end

network.Identifiers_kegg_reaction = network.reaction_KEGGID;


% -----------------------------------------------------------
% convert from M to mM

conc_min  = 1000 * cell_string2num(network.ConcentrationMin);
conc_max  = 1000 * cell_string2num(network.ConcentrationMax);

ind_water = find(strcmp('C00001',network.metabolite_KEGGID));
conc_min(ind_water) = 1;
conc_max(ind_water) = 1;

network                       = netgraph_make_graph(network);
network.graphics_par.metnames = network.metabolite_names;

network.reaction_formula_KEGGID = strrep(network.formulae,' ','_');

network = netgraph_read_positions(network,filenames.position_file,[0,0],1,1,network.reaction_formula_KEGGID);


% -----------------------------------------------------------
% make data structure 'kinetic_data' for all kinetic constants
% and insert values of equilibrium constants

import_quantity_list = {'standard chemical potential','Michaelis constant',...
                    'activation constant', 'inhibitory constant',...
                    'equilibrium constant','substrate catalytic rate constant', ...
                    'product catalytic rate constant'};

kinetic_data = data_integration_load_kinetic_data(import_quantity_list, [], network, kinetic_data_file_names, 0, 1,1,0,'Organism',organism_long);
