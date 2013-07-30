% --------------------------------------------------------------------------

obd_pathways_prepare_all;

% --------------------------------------------------------------------------

data_id  = 'd1_eco'; organism = 'eco'; obd_pathways_prepare_single;
data_id  = 'd2_bsu'; organism = 'bsu'; obd_pathways_prepare_single;
data_id  = 'd3_sce'; organism = 'sce'; obd_pathways_prepare_single;

% --------------------------------------------------------------------------

organism   = 'eco';  
model_name = 'EMP-GLYCOLYSIS-PTS'; 
run_id = 'r1_eco_KcatBrenda';   data_id  = 'd1_eco'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
run_id = 'r2_eco_NoKcat';       data_id  = 'd1_eco'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
run_id = 'r3_eco_KcatForward';  data_id  = 'd1_eco'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

organism   = 'bsu';  
model_name = 'EMP-GLYCOLYSIS-PTS'; 
run_id = 'r1_bsu_KcatBrenda';   data_id  = 'd2_bsu'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
run_id = 'r2_bsu_NoKcat';       data_id  = 'd2_bsu'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
run_id = 'r3_bsu_KcatForward';  data_id  = 'd2_bsu'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

organism   = 'sce';  
model_name = 'EMP-GLYCOLYSIS';
run_id = 'r1_sce_KcatBrenda';   data_id  = 'd3_sce'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
run_id = 'r2_sce_NoKcat';       data_id  = 'd3_sce'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
run_id = 'r3_sce_KcatForward';  data_id  = 'd3_sce'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

% --------------------------------------------------------------------------

data_id   = 'd1_eco'; run_id = 'FSC_eco_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;
data_id   = 'd2_bsu'; run_id = 'FSC_bsu_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;
data_id   = 'd3_sce'; run_id = 'FSC_sce_KcatBrenda'; kcat_usage = 'use'; script_fsc_obd_pathways;

script_fsc_fermentation_respiration
