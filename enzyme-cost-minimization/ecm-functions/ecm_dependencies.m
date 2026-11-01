function ecm_dependencies()

% ECM_DEPENDENCIES - Check whether the packages used by Enzyme Cost Minimization are installed
%
% ecm_dependencies()
%
% Called once per Matlab session by the main functions (see README.md for details)

if ~exist('matlab_utils_installed','file'),
  error('Please install the Matlab utility functions (https://github.com/liebermeister/matlab-utils)');
end

if ~exist('mnt_version','file'),
  error('Please install the Metabolic Network Toolbox (https://github.com/liebermeister/metabolic-network-toolbox)');
end

if exist('package_version_check','file')
  package_version_check('The Matlab utility functions', 'matlab_utils_version', '1.0.0');
  package_version_check('The Metabolic Network Toolbox', 'mnt_version', '1.3.0');
end

if ~exist('sbtab_toolbox_version','file'),
  warning('Please install SBtab for Matlab (https://github.com/liebermeister/sbtab-matlab) - Otherwise SBtab files cannot be read or written');
elseif exist('package_version_check','file')
  package_version_check('SBtab for Matlab', 'sbtab_toolbox_version', '1.0.0');
end

if ~exist('parameter_balancing_options','file'),
  warning('Please install Parameter Balancing (https://github.com/liebermeister/parameter-balancing-matlab) - It is required by ecm_parameter_balancing');
end

if ~exist('TranslateSBML','file'),
  warning('Please install the SBML Toolbox (http://sbml.org/Software/SBMLToolbox) - Otherwise SBML import and export will not work');
end
