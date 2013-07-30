% ------------------------------------------
% Sampling of metabolite concentrations compatible with predefined flux directions
% The linear constraint problem is explicitly formulated and solved by convex sampling 
%
% Input information: network structure; flux distribution; standard chemical potentials, 
% upper and lower bounds on concentrations; priors for individual metabolites 
%
% Note that concentrations are measured in Molar  (not in mM like in concentration_sampling_pb.m)!!
% Upper bound on A is not supported; for sampling, setting n_warm = n_sample
% ------------------------------------------

if 0,
 % filenames = sbtab_model_filenames('ecoli_glycolysis');       condition = []; concentration_sampling;
 % filenames = sbtab_model_filenames('ecoli_glycolysis2');      condition = []; concentration_sampling;
  filenames = sbtab_model_filenames('ecoli_entner_doudoroff'); condition = []; concentration_sampling;
  %%filenames = sbtab_model_filenames('ecoli_entner_doudoroff_pts'); condition = []; concentration_sampling;
  %%filenames = sbtab_model_filenames('bscm'); condition = 'glucose'; concentration_sampling;
  %%filenames = sbtab_model_filenames('ycm'); condition = 'diauxic_glucose'; concentration_sampling;
end

if 0,
 display('Using old results');
else,
  
options.conc_min    = 0.000001; % in molar
options.conc_max    = 0.1;      % in molar
options.flag_fix_concentrations = 1;
options.variability = 5;        % variability of known concentrations
options.n_sample    = 200; % 200
options.n_pick      = 50;  % 50
options.flag_artificial_centering = 0;
options.epsilon     = 0.5; % minimal Delta G (in kJ/mol) needed to drive a forward reaction


% ------------------------------------------
% load network model

load(filenames.network_file);  % network, network_CoHid 
load(filenames.flux_file);     % v

[nm,nr]   = size(network.N);


% ------------------------------------------
% load kinetic data directly from "Quantity" table (convention M)

kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','concentration'}, [], network, filenames.table_quantities, 0, 1);

mu0 = kinetic_data.mu0.mean; 

if isnan(mu0),
  error('mu0 values missing');
  mu0(isnan(mu0)) = 0;
end

% for which reactions are all chemical potentials known?
ind_ok     = find(abs(network.N') * isnan( kinetic_data.mu0.mean)==0);
N_sel      = network.N(:,ind_ok);
Keq_sel    = exp(- N_sel' * mu0/RT);
v_sign_sel = sign(v(ind_ok));


% ------------------------------------------

ind_irrel = label_names({'H+','H2O','Biomass'},network.metabolites);
ind_sel   = setdiff(1:nm,ind_irrel);
n_sel     = length(ind_sel);

log_c_min = log(options.conc_min) * ones(n_sel,1);
log_c_max = log(options.conc_max) * ones(n_sel,1);

% tighter constraints for known concentrations
ind_conc_known             = find(isfinite(kinetic_data.c.mean(ind_sel)));
log_c_min(ind_conc_known)  = log(kinetic_data.c.mean(ind_sel(ind_conc_known))) - log(options.variability);
log_c_max(ind_conc_known)  = log(kinetic_data.c.mean(ind_sel(ind_conc_known))) + log(options.variability);

% requirement: all fluxes in forward direction!
M_ineq = diag(v_sign_sel) * N_sel(ind_sel,:)';
b_ineq = diag(v_sign_sel) * log(Keq_sel) - options.epsilon/RT;

% fixed concentrations?
if options.flag_fix_concentrations, 
  ee   = eye(length(ind_sel));
  M_eq = ee(ind_conc_known,:);
  b_eq = log(kinetic_data.c.mean(ind_sel(ind_conc_known)));
else,
  M_eq = [];
  b_eq = [];
end


% -------------------------------------
% define general gaussian

n_sigmas   = 2;
log_c_mean = 1/2 * [log_c_max+log_c_min]; 
log_c_std  = 1/n_sigmas * 1/2 * [log_c_max-log_c_min];

% specific priors for known concentrations
ind_conc_known             = find(isfinite(kinetic_data.c.mean(ind_sel)));
log_c_mean(ind_conc_known) = log(kinetic_data.c.mean(ind_sel(ind_conc_known)));
log_c_std(ind_conc_known)  = log(options.variability);

log_c_cov  = diag(log_c_std.^2 );


% -------------------------------------
% run sampling

% uniform sampling with constraints
%X_sample_con = convex_sampling_hit_and_run_eqconst(M_ineq,b_ineq,M_eq,b_eq,log_c_min,log_c_max,[],[],options.n_sample,options.n_sample,options.n_pick,options.flag_artificial_centering);
X_sample_con = convex_sampling(M_ineq,b_ineq,M_eq,b_eq,options.n_sample,log_c_min,log_c_max);

% test: no flux constraints (good to see if enough samples are taken)
%X_sample_test = convex_sampling_hit_and_run_eqconst([],[],M_eq,b_eq,log_c_min,log_c_max,[],[],options.n_sample,options.n_sample,options.n_pick,options.flag_artificial_centering);
X_sample_test = convex_sampling([],[],M_eq,b_eq,options.n_sample,log_c_min,log_c_max);

% sampling with prior
X_sample_gau = convex_sampling_hit_and_run_eqconst(M_ineq,b_ineq,M_eq,b_eq,log_c_min,log_c_max,log_c_mean,log_c_cov,options.n_sample,options.n_sample,options.n_pick,options.flag_artificial_centering);

cd(filenames.model_dir); cd results
save('concentration_sampling_result','X_sample_con','X_sample_test','X_sample_gau','options','log_c_min','log_c_max','kinetic_data','ind_sel','network','network_CoHid','log_c_mean','ind_conc_known','v');

end

% -------------------------------------
% graphics

cd(filenames.model_dir); cd results
load('concentration_sampling_result');

mu0 = kinetic_data.mu0.mean; 

goptions.fig_offset = 0;
goptions.model_name = filenames.model_name;
goptions.conc_min = options.conc_min;
goptions.conc_max = options.conc_max;
goptions.network_CoHid = network_CoHid;
goptions.network_CoHid.graphics_par.squaresize = 0.03;
goptions.ind_conc_known  = ind_conc_known;
goptions.log_c_mean      = log_c_mean    ;
goptions.log_c_lower     = log_c_min     ;
goptions.log_c_upper     = log_c_max     ;
if options.flag_fix_concentrations,
  goptions.log_c_lower     = log_c_mean     ;
  goptions.log_c_upper     = log_c_mean     ;
end
goptions.v               = v     ;

goptions.run = 'sampling_constraint';
concentration_sampling_figures(X_sample_con,kinetic_data.mu0.mean,network,kinetic_data,ind_sel,goptions)

goptions.run = 'sampling_test';
concentration_sampling_figures(X_sample_test,kinetic_data.mu0.mean,network,kinetic_data,ind_sel,goptions)

goptions.run = 'sampling_gaussian';
concentration_sampling_figures(X_sample_gau,kinetic_data.mu0.mean,network,kinetic_data,ind_sel,goptions)


% -------------------------------------
% make and save network graphics

gpars = struct('actstyle','none','arrowstyle','fluxes','arrowsize',0.03,'actprintnames',1);
figure(100); clf; netgraph_concentrations(network,[],v,1,gpars)

cd  /home/wolfram/projekte/pathway_modelling/ps-files/concentration_sampling
print([ filenames.model_name '/' filenames.model_name '_network.eps'],'-f100','-depsc');


% ----------------------------------------------
% print data tables

% mu0 values (corresponding to concentrations in M)
print_matrix(kinetic_data.mu0.mean, network.metabolites)

% Delta mu0 values in direction of flux  (corresponding to concentrations in M)
print_matrix(sign(v) .* [network.N' * mu0], network.actions)

% concentrations (in M)
cc = kinetic_data.c.mean; 
print_matrix(cc(isfinite(cc)), network.metabolites(isfinite(cc)))
cd([filenames.model_dir '/results']);
[nm,nr] = size(network.N);
mu0 = kinetic_data.mu0.mean; 
dmu0 = sign(v) .* [network.N' * mu0];
ind_c = find(isfinite(kinetic_data.c.mean));
table([{'\textbf{Compound}','\textbf{mu0}','\textbf{Unit}'}; [network.metabolites, num2cell(kinetic_data.mu0.mean), repmat({'kJ/mol'},nm,1)]],'tex','concentration_sampling_data_mu0.tex')
table([{'\textbf{Compound}','\textbf{c}','\textbf{Unit}'}; [network.metabolites(ind_c), num2cell(kinetic_data.c.mean(ind_c)), repmat({'M'},length(ind_c),1)]],'tex','concentration_sampling_data_c.tex')
table([{'\textbf{Enzyme}','\textbf{dmu0}','\textbf{Unit}'}; [network.actions, num2cell(dmu0), repmat({'kJ/mol'},nr,1)]],'tex','concentration_sampling_data_dmu0.tex')
