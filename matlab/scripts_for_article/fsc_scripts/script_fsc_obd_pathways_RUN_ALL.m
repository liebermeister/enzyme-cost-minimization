% --------------------------------------------------------------------------
% Model set 'obd_pathways'
% --------------------------------------------------------------------------

% --------------------------------------------------------------------------
% prepare general files for the whole pathway set

pathway_collection = 'obd_pathways'; organisms = {'eco', 'bsu', 'sce'}; flag_thermo = 0; fsc_workflow_prepare_all;


% --------------------------------------------------------------------------
% prepare model-specific files

  data_id  = 'd1_eco'; pathway_collection = 'obd_pathways'; organism = 'eco'; fsc_workflow_prepare_single;
  data_id  = 'd2_bsu'; pathway_collection = 'obd_pathways'; organism = 'bsu'; fsc_workflow_prepare_single;
  data_id  = 'd3_sce'; pathway_collection = 'obd_pathways'; organism = 'sce'; fsc_workflow_prepare_single;
  data_id  = 'd0_eco'; pathway_collection = 'kegg_modules'; organism = 'eco'; fsc_workflow_prepare_single;


% --------------------------------------------------------------------------
% Script for a single pathway

% organism   = 'eco';  
% model_name = 'EMP-GLYCOLYSIS-PTS'; 
% run_id = 'r1_eco_KcatBrenda';   data_id  = 'd1_eco'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r2_eco_NoKcat';       data_id  = 'd1_eco'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r3_eco_KcatForward';  data_id  = 'd1_eco'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;
% 
% organism   = 'bsu';  
% model_name = 'EMP-GLYCOLYSIS-PTS'; 
% run_id = 'r1_bsu_KcatBrenda';   data_id  = 'd2_bsu'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r2_bsu_NoKcat';       data_id  = 'd2_bsu'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r3_bsu_KcatForward';  data_id  = 'd2_bsu'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;
% 
% organism   = 'sce';  
% model_name = 'EMP-GLYCOLYSIS';
% run_id = 'r1_sce_KcatBrenda';   data_id  = 'd3_sce'; kcat_usage = 'use';     script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r2_sce_NoKcat';       data_id  = 'd3_sce'; kcat_usage = 'none';    script_fsc_EMP_GLYCOLYSIS;
% run_id = 'r3_sce_KcatForward';  data_id  = 'd3_sce'; kcat_usage = 'forward'; script_fsc_EMP_GLYCOLYSIS;

% --------------------------------------------------------------------------
% Script for several pathways

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
data_id = 'd1_eco'; run_id = 'FSC_eco_KcatBrenda'; kcat_usage = 'use'; organism='eco';
model_names = {'EMP-GLYCOLYSIS-PTS-oxPPP-TCA'};
script_fsc_obd_pathways;
%%'TCA' 'ED-GLYCOLYSIS', 'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA'

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
data_id   = 'd2_bsu'; run_id = 'FSC_bsu_KcatBrenda'; kcat_usage = 'use'; organism='bsu';
model_names = {'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA'}; % 'TCA', 
script_fsc_obd_pathways;

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
data_id   = 'd3_sce'; run_id = 'FSC_sce_KcatBrenda'; kcat_usage = 'use'; organism='sce';
model_names = {'EMP-GLYCOLYSIS-TCA'}; %  'TCA', 'EMP-GLYCOLYSIS', 
script_fsc_obd_pathways;


% --------------------------------------------------------------------------
% Script for screening of flux modes

% script_fsc_fermentation_respiration
