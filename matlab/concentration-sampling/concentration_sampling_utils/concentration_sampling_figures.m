function concentration_sampling_figures(log_c_matrix, mu0_matrix, network, kinetic_data, ind_rel, options)

default_options = struct('print_flag','1');
options = join_struct(default_options, options);


% log_c_matrix,mu0_matrix refer only to the relevant metabolites

if ~isfield(options,'network_CoHid'), options,network_CoHid = network; end; 

if size(mu0_matrix,2) ==1, 
  mu0_matrix = repmat(  mu0_matrix,1,size(log_c_matrix,2));
end

[nm,nr] = size(network.N);
[n_rel,n_sample] = size(log_c_matrix);

mu_matrix            =   mu0_matrix;
mu_matrix(ind_rel,:) =   mu0_matrix(ind_rel,:) + RT * log_c_matrix;
A_matrix  = - network.N' * mu_matrix;

figure(options.fig_offset+1); clf 
cs_display_cumulative_distribution(options.fig_offset+1,exp(log_c_matrix),standardise_metabolite_names(network.metabolites(ind_rel)),'Concentrations (cumulative distribution)',options.conc_min, options.conc_max,struct('flag_logarithmic',1));

my_A_matrix = diag(sign(options.v)) * A_matrix;
figure(options.fig_offset+2); clf 
cs_display_cumulative_distribution(options.fig_offset+2,my_A_matrix,network.actions,'Reaction affinities along flux (cumulative distribution)',0,60,options);

figure(options.fig_offset+3); clf 
cs_display_distribution_errorbar(options.fig_offset+3,exp(log_c_matrix),standardise_metabolite_names(network.metabolites(ind_rel)),'Concentrations (median and 5% quantiles)',options.conc_min, options.conc_max, kinetic_data.c.median(ind_rel),struct('flag_logarithmic',1));

log_c_median = median(log_c_matrix')';
figure(options.fig_offset+4); clf; 
log_c_all_median = nan*zeros(nm,1);
log_c_all_median(ind_rel) = log_c_median;
netgraph_concentrations(options.network_CoHid,log_c_all_median,[],1,struct('actprintnames',0,'metvaluesmin',log(options.conc_min),'metvaluesmax',log(options.conc_max),'showsign',0));
title('Median log concentrations (color scale spans allowed range)');


figure(options.fig_offset+5);
mu_median          = median(mu_matrix')'; 
A_median           = median(A_matrix')';
netgraph_concentrations(options.network_CoHid,mu_median,[],1,struct('actprintnames',0,'arrowvalues',A_median,'arrowstyle','fluxes','arrowsize',0.05));
title('Median chemical potentials and reaction affinities');


figure(options.fig_offset+6);
mu0_median =   median(mu0_matrix')';
A0_matrix  = - network.N' * mu0_matrix;
A0_median  = median(A0_matrix')';
netgraph_concentrations(options.network_CoHid,mu0_median,A0_median.*sign(options.v),1,struct('actprintnames',1,'arrowvalues',A0_median,'arrowstyle','fluxes','arrowsize',0.05,'actprintvalues',1,'metprintvalues',0,'metstyle','none','actvaluesmin',-30,'actvaluesmax',30));
title('Median std. chemical potentials and standard reaction affinities');


if length( options.print_flag),
  cd(['/home/wolfram/projekte/pathway_modelling/ps-files/concentration_sampling/' options.model_name ]);
  print([  options.model_name '_' options.run '_cum.eps'],    ['-f' num2str(options.fig_offset+1)],'-depsc');
  print([  options.model_name '_' options.run '_aff_cum.eps'],['-f' num2str(options.fig_offset+2)],'-depsc');
  print([  options.model_name '_' options.run '_bar.eps'],    ['-f' num2str(options.fig_offset+3)],'-depsc');
  print([  options.model_name '_' options.run '_net.eps'],    ['-f' num2str(options.fig_offset+4)],'-depsc');
  print([  options.model_name '_' options.run '_aff_net.eps'],['-f' num2str(options.fig_offset+5)],'-depsc');
  print([  options.model_name '_' options.run '_mu0_net.eps'],['-f' num2str(options.fig_offset+6)],'-depsc');
end
