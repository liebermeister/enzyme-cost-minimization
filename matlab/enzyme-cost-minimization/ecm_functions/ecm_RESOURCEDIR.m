% ECM_RESOURCEDIR - Return ECM resource directory
%
% d = ecm_RESOURCEDIR()

function d = ecm_RESOURCEDIR()

d = [ecm_BASEDIR '..' filesep '..' filesep 'resources' filesep];
