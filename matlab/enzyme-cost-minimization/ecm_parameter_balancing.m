function [r, r_orig, kinetic_data, r_samples, ecm_options, parameter_prior, r_std] = ecm_parameter_balancing(network, ecm_options, kinetic_data);

% ECM_PARAMETER_BALANCING - Prepare and run parameter balancing with standard options
%
% [r, r_orig, kinetic_data, r_samples, ecm_options, parameter_prior, r_std] = ecm_parameter_balancing(network, ecm_options, kinetic_data);
%
% Output
%   r        Kinetic constants (posterior model values; used as input in parameter balancing)
%   r_orig   Original kinetic constants (used as input in parameter balancing)
%   r_std    Kinetic constants (posterior standard deviations)
% 
% Uses (potentially) the following options from options.
%   options.flag_given_kinetics
%   options.reaction_column_name  (only if no kinetic data are given)
%   options.compound_column_name  (only if no kinetic data are given)
%   options.kcat_usage            {'use','none','forward'} (default: 'use')
%   options.kcat_prior_median
%   options.kcat_prior_log10_std
%   options.kcat_lower
%   options.kcatr_lower
%   options.kcat_upper
%   options.KM_lower
%   options.Keq_upper
%   options.parameter_prior_file
%   options.GFE_fixed
%   options.use_pseudo_values
%   options.fix_Keq_in_sampling

ecm_options = join_struct(ecm_default_options(network),ecm_options);

% legacy: this has been used in EFM paper and elsewhere!
ecm_options.use_pseudo_values = 0; 


% --------------------------------------------------------------------------------------

if ecm_options.flag_given_kinetics,
  
  %% If desired, use kinetic parameters directly from the model
  switch network.kinetics.type 
    case {'cs','ms','ds'},
      display('  Using kinetic constants found in network.kinetics');
      r.Keq     = network.kinetics.Keq;
      r.KM      = network.kinetics.KM;
      [r.Kcatf, r.Kcatr] = modular_KV_Keq_to_kcat(network.N,network.kinetics);
    otherwise error('Kinetics cannot be handled');
  end

else
  
  %% Use kinetic parameters from "kinetic_data" function argument and run parameter balancing
  [r, r_orig, kinetic_data, r_samples, parameter_prior, r_std] = parameter_balancing_kinetic(network, kinetic_data,ecm_options);

end

% -----------------------------
% if desired, insert original Keq values

if ecm_options.insert_Keq_from_data,
  if length(Keq_given),
    display('Using predefined equilibrium constants exactly'); 
    ind_finite = find(isfinite(Keq_given));
    r.Keq(ind_finite) = r_orig.Keq(ind_finite);
  end
end
