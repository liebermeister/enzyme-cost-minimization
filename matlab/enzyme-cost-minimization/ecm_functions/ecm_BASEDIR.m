% ECM_BASEDIR - Return ECM toolbox directory
%
% d = ecm_BASEDIR()


function d = ecm_BASEDIR()

d = [fileparts(which(mfilename)) filesep '..' filesep];
