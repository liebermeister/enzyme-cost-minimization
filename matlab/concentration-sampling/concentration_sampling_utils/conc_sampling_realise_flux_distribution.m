function [c,Keq,u,network,A] = conc_sampling_realise_flux_distribution(network,c,v,u)

A_min = 0.1; % minimal force to drive a reaction

network.kinetics.c = c;
network.kinetics.u = u;
modular_rate_law_to_sbtab(network,'/tmp/dummi.tsv',struct('use_sbml_ids',0)); % export values
data_quantities      = {'concentration','equilibrium constant','reaction affinity'}';
%kinetic_data         = data_integration_load_kinetic_data(data_quantities, [], network, {'/tmp/dummi.tsv', '/home/wolfram/projekte/pathway_modelling/models/bsu_central_metabolism/bsu_central_metabolism_QuantityData.tsv'}, 0, 0);
kinetic_data         = data_integration_load_kinetic_data(data_quantities, [], network, {'/tmp/dummi.tsv'}, 0, 0);
kinetic_data.Keq.std = 0.0001*kinetic_data.Keq.std;
kinetic_data.c.std(find(network.external)) = 0.01*kinetic_data.c.std(find(network.external));
kinetic_data.A.lower(v>0) = A_min;
kinetic_data.A.upper(v<0) = -A_min;
model_quantities  = {'equilibrium constant','concentration','reaction affinity'}';
basic_quantities  = {'standard chemical potential','concentration'}';
task              = parameter_balancing_task(network, kinetic_data, [], model_quantities, basic_quantities);
res               = parameter_balancing(task, [], struct('insert_pseudo_values',0));

c   = res.kinetics_posterior_mode.c;
Keq = res.kinetics_posterior_mode.Keq;
A = res.kinetics_posterior_mode.A;
network.kinetics.c   = c;
network.kinetics.Keq = Keq;
my_v = network_velocities(c,network);
u = v./my_v;

network.kinetics.u = u;
