function setup_test_path()

% setup_test_path()
%
% Add the package to the Matlab path and check that the required packages are installed
% (called at the start of each test file)

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'enzyme-cost-minimization')));

ecm_dependencies;
