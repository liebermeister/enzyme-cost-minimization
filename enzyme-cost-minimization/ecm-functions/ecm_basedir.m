function d = ecm_basedir()

% ECM_BASEDIR - Return the main directory of the Enzyme Cost Minimization repository
%
% d = ecm_basedir()

d = fileparts(fileparts(fileparts(which(mfilename))));
