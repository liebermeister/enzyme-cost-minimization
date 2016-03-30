function [network, network_split, network_simple] = gams_infiles_to_network(directory) 

% directory = '/home/wolfram/projekte/flux_cost_functions/gams/resultsCarlson/network3-50';
% [network, network_split, network_simple] = gams_infiles_to_network(directory) 
% save meike_network network network_split network_simple

cd(directory)
a = load_any_table('kmstoich.csv',',');

metabolites = unique(a(:,1));
reactions   = unique(a(:,2));

N = zeros(length(metabolites), length(reactions));

i1 = label_names(a(:,1),metabolites);
i2 = label_names(a(:,2),reactions);
i3 = cell_string2num(a(:,3));

for it = 1:length(i1);
  N(i1(it),i2(it))= i3(it);
end

b = load_any_table('metfixed.csv',',');
ll = label_names(b(:,1),metabolites);
ll = ll(find(ll));

network = network_construct(N,zeros(size(reactions)),ll,metabolites,reactions);
network = netgraph_make_graph(network);

network_split  = netgraph_split_metabolites(network,{'atp','adp','nad','nadh','co2'});
network_simple = netgraph_simple_graph(network,{'atp','adp','nad','nadh','co2'});

