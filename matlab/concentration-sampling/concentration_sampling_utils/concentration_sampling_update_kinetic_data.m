function kinetic_data = concentration_sampling_update_kinetic_data(kinetic_data, network, v, options)

% kinetic_data = concentration_sampling_update_kinetic_data(kinetic_data, network, v, options)
%
% insert some data into 'kinetic_data' data structure and make sure it's consistent

options_default.variability = 0.1;
options_default.A_mean      = 5;
options_default.A_std       = 0.5;

options = join_struct(options_default,options);


[nm,nr] = size(network.N);

% ------------------------------------------
% set bounds on metabolite concentrations

kinetic_data.c.lower  = options.conc_min * ones(nm,1);
kinetic_data.c.upper  = options.conc_max * ones(nm,1);


% ------------------------------------------
% insert predefined variability for known metabolites

if isfinite(options.variability),
  ind = find(isfinite(kinetic_data.c.median));
  kinetic_data.c.lower(ind)    = kinetic_data.c.median(ind) / options.variability;
  kinetic_data.c.upper(ind)    = kinetic_data.c.median(ind) * options.variability;
  kinetic_data.c.lower_ln(ind) = log(kinetic_data.c.lower(ind));
  kinetic_data.c.upper_ln(ind) = log(kinetic_data.c.upper(ind));
  kinetic_data.c.mean_ln(ind)  = 0.5 * [kinetic_data.c.upper_ln(ind) + kinetic_data.c.lower_ln(ind)]; 
  kinetic_data.c.std_ln(ind)   = 0.5 * [kinetic_data.c.upper_ln(ind) - kinetic_data.c.lower_ln(ind)]; 
  [kinetic_data.c.mean(ind), kinetic_data.c.std(ind)] = lognormal_log2normal(kinetic_data.c.mean_ln(ind), kinetic_data.c.std_ln(ind));
end

% ------------------------------------------
% fix water at virtual standard concentration 1mM

ind_water = label_names({'H2O'},network.metabolites);
if ind_water,
  kinetic_data.c.mean(ind_water)  = 1;
  kinetic_data.c.std(ind_water)   = 0.001;
  kinetic_data.c.lower(ind_water) = 0.99;
  kinetic_data.c.upper(ind_water) = 1.01;
end

% ------------------------------------------
% make values consistent

kinetic_data.c.lower_ln = log(kinetic_data.c.lower);
kinetic_data.c.upper_ln = log(kinetic_data.c.upper);
[kinetic_data.c.mean_ln, kinetic_data.c.std_ln] = lognormal_normal2log(kinetic_data.c.mean, kinetic_data.c.std);  
kinetic_data.c.median = exp(kinetic_data.c.mean_ln);


% ------------------------------------------
% insert missing mu0 standard deviations

kinetic_data.mu0.std   = options.sigma_mu0 * ones(size(kinetic_data.mu0.std));
kinetic_data.mu0.lower = kinetic_data.mu0.median - options.sigma_mu0;
kinetic_data.mu0.upper = kinetic_data.mu0.median + options.sigma_mu0;

% ------------------------------------------
% set bounds on reaction affinities

ind_A = find(~isfinite(kinetic_data.A.mean));
kinetic_data.A.mean(ind_A) = options.A_mean * sign(v(ind_A));
kinetic_data.A.std(ind_A)  = options.A_std  * ones(length(ind_A),1);
kinetic_data.A.median      = kinetic_data.A.mean;

ind_A = find(~isfinite(kinetic_data.A.lower));
kinetic_data.A.lower(ind_A) = - options.A_max * ones(length(ind_A),1);
kinetic_data.A.upper(ind_A) =   options.A_max * ones(length(ind_A),1);

kinetic_data.A.lower(v>0) =  options.A_min;
kinetic_data.A.upper(v<0) = -options.A_min;

% ------------------------------------------
% set range for enzyme levels

if isfield(options,'u_min'),
  kinetic_data.u.lower    = options.u_min * ones(nr,1);
  kinetic_data.u.upper    = options.u_max * ones(nr,1);
  kinetic_data.u.lower_ln = log(kinetic_data.u.lower);
  kinetic_data.u.upper_ln = log(kinetic_data.u.upper);
end