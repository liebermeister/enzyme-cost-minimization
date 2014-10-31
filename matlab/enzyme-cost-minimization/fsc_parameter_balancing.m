function [r, r_orig, kinetic_data] = fsc_parameter_balancing(network,fsc_options, kinetic_data);

% [r, r_orig, kinetic_data] = fsc_parameter_balancing(network,fsc_options, kinetic_data);
%
% FSC_PARAMETER_BALANCING - Prepare and run parameter balancing
%
% Output
%   r             Kinetic constants (used as input in parameter balancing)
%   r_orig        Original kinetic constants (used as input in parameter balancing)
% 
% Uses (potentially) the following options from fsc_options.
%  fsc_options.flag_given_kinetics
%  fsc_options.reaction_column_name
%  fsc_options.compound_column_name
%  fsc_options.kcat_usage
%  fsc_options.kcat_prior_median
%  fsc_options.kcat_prior_log10_std
%  fsc_options.kcat_lower
%  fsc_options.kcat_upper
%  fsc_options.KM_lower
%  fsc_options.Keq_upper
%  fsc_options.quantity_info_file


fsc_options_default = fsc_default_options(network);
fsc_options         = join_struct(fsc_options_default,fsc_options);


% -------------------------------------------------
% fsc_options

if fsc_options.flag_given_kinetics == 0,
  if isempty(kinetic_data),
    kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','standard chemical potential difference','Michaelis constant','activation constant',  'inhibitory constant','equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}, [], network, [], 0, 1, fsc_options.reaction_column_name, fsc_options.compound_column_name);
  end
end


% --------------------------------------------------------------------------------------
% run parameter balancing for kinetic constants

if fsc_options.flag_given_kinetics,

  switch network.kinetics.type 
    case {'cs','ms','ds'},
      display('Using kinetic constants found in network.kinetics');
      r.Keq     = network.kinetics.Keq;
      r.KM      = network.kinetics.KM;
      [r.Kcatf, r.Kcatr] = modular_KV_Keq_to_kcat(network.N,network.kinetics);
    otherwise error('Kinetics cannot be handled');
  end
  
else
  
  switch fsc_options.kcat_usage
    case 'use',
    
    case 'none',

      %% do not use any kcat data
      kk = kinetic_data;
      em = nan * kk.KM.median;
      kk.KM.median = em; kk.KM.mean = em; kk.KM.std = em; kk.KM.mean_ln = em; kk.KM.std_ln = em;
      em = nan * kk.KM.median;
      kk.KA.median = em; kk.KA.mean = em; kk.KA.std = em; kk.KA.mean_ln = em; kk.KA.std_ln = em;
      em = nan * kk.KI.median;
      kk.KI.median = em; kk.KI.mean = em; kk.KI.std = em; kk.KI.mean_ln = em; kk.KI.std_ln = em;
      em = nan * kk.Kcatf.median;
      kk.Kcatf.median = em; kk.Kcatf.mean = em; kk.Kcatf.std = em; kk.Kcatf.mean_ln= em; kk.Kcatf.std_ln = em;
      em = nan * kk.Kcatr.median;
      kk.Kcatr.median = em; kk.Kcatr.mean = em; kk.Kcatr.std = em; kk.Kcatr.mean_ln= em; kk.Kcatr.std_ln = em;
      kinetic_data = kk;
      
    case 'forward',
      %% invent kcat data such that forward values (along flux) ~= 
      %% match the standard value
      ind_p = find([v>=0]+isnan(v));
      ind_m = find(v<0);
      emp   = ones(size(ind_p));
      emm   = ones(size(ind_m));
      
      if isempty(fsc_options.kcat_prior_median), error('Kcat standard value missing'); end
      kcat_forward_value      = fsc_options.kcat_prior_median; % unit: 1/s
      kk                      = kinetic_data;
      kk.Kcatf.median(ind_p)  =     kcat_forward_value  * emp; 
      kk.Kcatr.median(ind_m)  =     kcat_forward_value  * emm; 
      kk.Kcatf.mean_ln(ind_p) = log(kcat_forward_value) * emp; 
      kk.Kcatr.mean_ln(ind_m) = log(kcat_forward_value) * emm; 
      kk.Kcatf.std_ln(ind_p)  = 0.1 * emp;
      kk.Kcatr.std_ln(ind_m)  = 0.1 * emm;
      [kk.Kcatf.mean(ind_p),kk.Kcatf.std(ind_p)] = lognormal_log2normal(kk.Kcatf.mean_ln(ind_p),kk.Kcatf.std_ln(ind_p));
      [kk.Kcatr.mean(ind_m),kk.Kcatr.std(ind_m)] = lognormal_log2normal(kk.Kcatr.mean_ln(ind_m),kk.Kcatr.std_ln(ind_m));
      kk.Kcatf.mean(ind_m) = nan;
      kk.Kcatr.mean(ind_p) = nan;
      kk.Kcatf.std(ind_m)  = nan;
      kk.Kcatr.std(ind_p)  = nan;
      kinetic_data         = kk;
  end
    
  [r, r_orig, kinetic_data] = parameter_balancing_kinetic(network, kinetic_data,[],[],struct('kcat_prior_median',fsc_options.kcat_prior_median,'kcat_prior_log10_std',fsc_options.kcat_prior_log10_std,'kcat_lower',fsc_options.kcat_lower,'kcat_upper',fsc_options.kcat_upper,'KM_lower',fsc_options.KM_lower,'Keq_upper',fsc_options.Keq_upper,'GFE_fixed',1,'quantity_info_file',fsc_options.quantity_info_file));

end
