function d = ecm_resourcedir()

% ECM_RESOURCEDIR - Return the directory 'resources' (example models and data) of the repository
%
% d = ecm_resourcedir()

d = [ecm_basedir filesep 'resources'];
