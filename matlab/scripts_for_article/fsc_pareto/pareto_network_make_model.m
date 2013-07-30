% -----------------------------------------
% make example model

if 0,
  %% small branched network

  N = [-1  0  0  0; ...  
       1 -1 -1 -1; ...  
       0  1  0  0; ...
       0  0  1  0; ...
       0  0  0  1];

  ind_ext = [1 3 4 5]';

  network = network_construct(N,ones(4,1), ind_ext);
  network.kinetics = set_kinetics(network,'cs');
  network = netgraph_make_graph(network);
  network = netgraph_move(network);   
  
  cd ~/projekte/pathway_modelling/flux_specific_cost/matlab_flux_specific_cost/pareto
  save pareto_branch_network network
end
