% ------------------------------------------
% Sampling of metabolite concentrations compatible with predefined flux directions
% The problem is formulated and solved as a parameter balancing problem 
%
% Input information: network structure; flux distribution; standard chemical potentials, 
% upper and lower bounds on concentrations; priors for individual metabolites 
%
% Note that concentrations are measured in mM (not in M like in concentration_sampling.m)!!
%
% ------------------------------------------

if 0,
  filenames = sbtab_model_filenames('ecoli_glycolysis');       condition = []; concentration_sampling_pb;
  filenames = sbtab_model_filenames('ecoli_glycolysis2');      condition = []; concentration_sampling_pb;
  filenames = sbtab_model_filenames('ecoli_entner_doudoroff'); condition = []; concentration_sampling_pb;
  filenames = sbtab_model_filenames('MOGtranshydrogenase1');   condition = []; concentration_sampling_pb;
  filenames = sbtab_model_filenames('bscm'); condition = 'glucose';         concentration_sampling_pb;
  filenames = sbtab_model_filenames('ycm');  condition = 'diauxic_glucose'; concentration_sampling_pb;
end


options.conc_min    = 0.001; % in mM
options.conc_max    = 100;   % in mM
options.variability = 2;     % variability of known concentrations
options.A_max       = 60;
options.A_min       = 0.5;   % minimal Delta G needed to drive a forward reaction
options.A_mean      = nan;   % mean and std values for forward directions
options.A_std       = nan;
options.sigma_mu0   = 3;     % error of mu0 values (kJ/mol); 3 for alberty data

options.n_sample    = 100;   % 200
options.n_pick      = 100;   % 200 
options.n_warm      = [];
options.q_width     = 20;
options.try_rejection_sampling = 1;

switch filenames.model_name, 
  case 'MOGtranshydrogenase1',
    options.try_rejection_sampling  = 0;
end

options.start_sampling_near_mode    = 1;
options.flag_artificial_centering   = 1;
options.compute_reaction_velocities = 0;
options.data_refer_to_molar         = 1;


% ------------------------------------------

load(filenames.network_file);  % network, network_CoHid
load(filenames.flux_file);     % v

[nm,nr]   = size(network.N);
network.kinetics = set_cs_kinetics(network);

ind_water = label_names({'H2O'},network.metabolites);


% ------------------------------------------
% load preprocessed kinetic data (convention: mM) from several general files

% data (convention: mM) constructed from several general files
% load(filenames.kinetic_quantities_mat); % kinetic_data 

% data from "Quantity" table needs to be converted to mM first

kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','concentration','reaction affinity'}, [], network, filenames.table_quantities, 0, 1);

if options.data_refer_to_molar,

  ind_cor = setdiff(1:nm,ind_water);

  kinetic_data.mu0.median(ind_cor)  = kinetic_data.mu0.median(ind_cor) - RT * log(1000); 
  kinetic_data.mu0.mean(ind_cor)    = kinetic_data.mu0.mean(ind_cor)   - RT * log(1000); 
  kinetic_data.mu0.lower(ind_cor)   = kinetic_data.mu0.lower(ind_cor)  - RT * log(1000); 
  kinetic_data.mu0.upper(ind_cor)   = kinetic_data.mu0.upper(ind_cor)  - RT * log(1000); 
  kinetic_data.c.median             = kinetic_data.c.median   * 1000; 
  kinetic_data.c.mean               = kinetic_data.c.mean     * 1000; 
  kinetic_data.c.std                = kinetic_data.c.std      * 1000; 
  kinetic_data.c.lower              = kinetic_data.c.lower    * 1000; 
  kinetic_data.c.upper              = kinetic_data.c.upper    * 1000; 
  kinetic_data.c.mean_ln            = kinetic_data.c.mean_ln  + log(1000); 
  kinetic_data.c.lower_ln           = kinetic_data.c.lower_ln + log(1000); 
  kinetic_data.c.upper_ln           = kinetic_data.c.upper_ln + log(1000); 
  
  kinetic_data = concentration_sampling_update_kinetic_data(kinetic_data, network, v, options);

end

% -----------------------
% fix mu0 values and concentrations exactly at data values:

kinetic_data.mu0.std(find(isfinite(kinetic_data.mu0.std)))   = 10^-5;
kinetic_data.c.std(find(isfinite(kinetic_data.c.std)))       = 10^-5;
kinetic_data.c.std_ln(find(isfinite(kinetic_data.c.std_ln))) = 10^-5;


% ------------------------------------------
% run parameter balancing

quantity_info     = data_integration_load_quantity_info;
model_quantities  = {'standard chemical potential','concentration','reaction affinity'}';
data_quantities   = {'standard chemical potential','concentration','reaction affinity'}';
basic_quantities  = {'standard chemical potential','concentration'}';

%model_quantities  = {'standard chemical potential','concentration','reaction affinity', 'substrate catalytic rate constant', 'product catalytic rate constant', 'Michaelis constant'}';
%data_quantities   = {'standard chemical potential','concentration','reaction affinity', 'substrate catalytic rate constant', 'product catalytic rate constant','Michaelis constant','equilibrium constant'}';
%basic_quantities  = {'standard chemical potential','concentration', 'catalytic rate constant geometric mean', 'Michaelis constant'}';

task = parameter_balancing_task(network, kinetic_data, quantity_info, model_quantities, basic_quantities);
res  = parameter_balancing(task, quantity_info, struct('insert_pseudo_values',0));


% ------------------------------------------
% sample from parameter balancing results

options.insert_pseudo_values = 0;
options.x_start              = res.q_posterior.mode;

%% uniform sampling

options.sampling_type = 'uniform';

[q_sample_uni, xmodel_sample_uni, xdata_sample_uni, kinetics_sample_uni, v_sample_uni, q_min_uni, q_max_uni] = parameter_balancing_constrained_sampling(task, res, options, quantity_info);

%% gaussian sampling

options.sampling_type = 'gaussian';

[q_sample_gau, xmodel_sample_gau, xdata_sample_gau, kinetics_sample_gau, v_sample_gau, q_min_gau, q_max_gau] = parameter_balancing_constrained_sampling(task, res, options, quantity_info);

%% test: uniform sampling without thermodynamic constraints

kk          = kinetic_data;
kk.A.lower  = nan * ones(nr,1);
kk.A.upper  = nan * ones(nr,1);
kk.A.mean   = nan * ones(nr,1);
kk.A.std    = nan * ones(nr,1);

qq          = quantity_info;
qq.LowerBound{qq.symbol_index.A} = '-inf';
qq.UpperBound{qq.symbol_index.A} = 'inf'; 

task                  = parameter_balancing_task(network, kk, qq, model_quantities, basic_quantities);
options.sampling_type = 'uniform';
[q_sample_test, xmodel_sample_test, xdata_sample_test, kinetics_sample_test, v_sample_test, q_min_test, q_max_test] = parameter_balancing_constrained_sampling(task, res, options, qq);


% ------------------------------------------
% save results

cd(filenames.model_dir); cd results

save('concentration_pb_sampling_result','q_sample_uni', 'xmodel_sample_uni', 'xdata_sample_uni', 'kinetics_sample_uni' ,'q_sample_gau', 'xmodel_sample_gau', 'xdata_sample_gau', 'kinetics_sample_gau', 'q_sample_test', 'xmodel_sample_test', 'xdata_sample_test', 'kinetics_sample_test' ,'task' , 'res', 'options', 'quantity_info','network','network_CoHid','kinetic_data','v');


% ------------------------------------------
% graphics for posterior mode

ca

cd(filenames.model_dir); cd results
load('concentration_pb_sampling_result');

[nm,nr] = size(network.N);

figure(1); clf; 
netgraph_concentrations(network_CoHid,log10(kinetic_data.c.median),[],1,struct('actprintnames',0));
title('Median concentrations (data)');

figure(2); clf; 
netgraph_concentrations(network_CoHid,log10(res.kinetics_posterior_mode.c),[],1,struct('actprintnames',0));
title('Concentrations (posterior mode)');

figure(3); clf; 
netgraph_concentrations(network_CoHid,[],res.kinetics_posterior_mode.A,1,struct('actstyle','none','arrowstyle','fluxes','arrowsize',0.04,'actprintnames',0));
title('Reaction affinities (posterior mode)');

cd /home/wolfram/projekte/pathway_modelling/ps-files/concentration_sampling
print([filenames.model_name '/' filenames.model_name '_sampling_pb_c_data_median.eps'], '-f1', '-depsc');
print([filenames.model_name '/' filenames.model_name '_sampling_pb_c_post_mode.eps'], '-f2', '-depsc');
print([filenames.model_name '/' filenames.model_name '_sampling_pb_A_post_mode.eps'], '-f3', '-depsc');


% ------------------------------------------
% graphics for sampling results

goptions.fig_offset    = 0;
goptions.model_name    = filenames.model_name;
goptions.conc_min      = options.conc_min;
goptions.conc_max      = options.conc_max;
goptions.network_CoHid = network_CoHid;
goptions.ind_conc_known= find(isfinite(kinetic_data.c.mean));
goptions.log_c_mean    = kinetic_data.c.mean_ln ;
goptions.log_c_lower   = kinetic_data.c.lower_ln;
goptions.log_c_upper   = kinetic_data.c.upper_ln;
goptions.v             = v;

goptions.run           = 'sampling_pb_uni';
concentration_sampling_figures(log(kinetics_sample_uni.c),kinetics_sample_uni.mu0,network,kinetic_data,1:nm,goptions)

goptions.run           = 'sampling_pb_gau';
concentration_sampling_figures(log(kinetics_sample_gau.c),kinetics_sample_gau.mu0,network,kinetic_data,1:nm,goptions)

goptions.run           = 'sampling_pb_test';
concentration_sampling_figures(log(kinetics_sample_test.c),kinetics_sample_test.mu0,network,kinetic_data,1:nm,goptions)


% -------------------------------------
% make and save network graphics

gpars = struct('actstyle','none','arrowstyle','fluxes','arrowsize',0.03,'actprintnames',1);

figure(100); clf; 
netgraph_concentrations(network,[],v,1,gpars)

cd  /home/wolfram/projekte/pathway_modelling/ps-files/concentration_sampling
print([ filenames.model_name '/' filenames.model_name '_network.eps'],'-f100','-depsc');


% ----------------------------------------------
% print data tables

cd([filenames.model_dir '/results']);

[nm,nr] = size(network.N);
mu0     = kinetic_data.mu0.mean; 
dmu0    = sign(v) .* [network.N' * kinetic_data.mu0.mean];
ind_c   = find(isfinite(kinetic_data.c.mean));

table([{'\textbf{Compound}','\textbf{mu0}','\textbf{Unit}'}; [network.metabolites, num2cell(kinetic_data.mu0.mean), repmat({'kJ/mol'},nm,1)]],'tex','concentration_sampling_pb_data_mu0.tex')

table([{'\textbf{Compound}','\textbf{c}','\textbf{Unit}'}; [network.metabolites(ind_c), num2cell(kinetic_data.c.mean(ind_c)),repmat({'mM'},length(ind_c),1)]],'tex','concentration_sampling_pb_data_c.tex')

table([{'\textbf{Enzyme}','\textbf{dmu0}','\textbf{Unit}'}; [network.actions, num2cell(dmu0), repmat({'kJ/mol'},nr,1)]],'tex','concentration_sampling_pb_data_dmu0.tex')
