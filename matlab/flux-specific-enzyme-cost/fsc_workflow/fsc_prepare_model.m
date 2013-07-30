function [network, filenames, options] = fsc_prepare_model(model_directory)

% [network, filenames, options] = fsc_prepare_model(model_directory)
%
% Prepare modelfor FSC optimisation: 
% Read model in SBtab format and create .mat and SBML files

% TEST
if 0,
  model_directory = '~/projekte/flux_specific_cost/fsc/models/ecoli_glycolysis/';
  [network, filenames, options] = fsc_prepare_model(model_directory);
end


% -----------------------------------------------
% read options file 'fsc_options.txt' from model directory
% and define filenames

default_options = fsc_file_locations; 
options         = fsc_read_options([ model_directory '/fsc_options.txt'], fsc_file_locations);
filenames       = fsc_filenames(model_directory, model_directory,options.model_name);


%--------------------------------------
% import network from SBtab files

model_sbtab = [model_directory '/' options.model_name '_network'];
opt         = struct('load_quantity_table',0);
network     = sbtab_to_network(model_sbtab, opt);


%--------------------------------------
% graphics

network = netgraph_make_graph(network);
network = netgraph_read_positions(network,filenames.table_positions,[0,0],1);

if 0,
  %% Move elements interactively
  network = netgraph_move(network,'text');
  netgraph_print_positions(network,filenames.table_positions);
end

network_CoHid = network;
netgraph_draw(network);


%--------------------------------------
% save network (as .mat file and SBML file)

save(filenames.network_file,'network','network_CoHid','options');

network2sbml(network,0,options.model_name,filenames.sbml_file);

%network_sbml_export(network,0,options.model_name,filenames.sbml_file,'',2,1)


%--------------------------------------
% compute fluxes (and save as .mat file)

K = full(network_analyse(network));
v = K(:,1) * sign(K(1));

%gpars = struct('actstyle','none','arrowstyle','fluxes','arrowsize',0.03,'actprintnames',0);
%figure(100); clf; netgraph_concentrations(network,[],sign(v),1,gpars)

save(filenames.flux_file, 'v');

display(sprintf('Writing files to directory "%s"',model_directory));