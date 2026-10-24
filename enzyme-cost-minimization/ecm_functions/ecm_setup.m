% ECM_SETUP - Determine data directory

function ecm_info = ecm_setup()

ecm_info.data_dir = [ecm_BASEDIR '..' filesep 'resources' filesep 'data-files' filesep];
