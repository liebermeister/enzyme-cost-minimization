% -------------------------------------------------------------------
% Display all results of script_fsc_obd_pathways in one go
% -------------------------------------------------------------------

if 0,
  data_id = 'd1_eco'; run_id = 'FSC_eco_KcatBrenda'; kcat_usage = 'use'; organism='eco';
  model_names = {'EMP-GLYCOLYSIS-PTS-oxPPP-TCA'};% , 'EMP-GLYCOLYSIS-PTS'};
  script_fsc_obd_pathways_display;

  data_id   = 'd2_bsu'; run_id = 'FSC_bsu_KcatBrenda'; kcat_usage = 'use'; organism='bsu';
  model_names = {'EMP-GLYCOLYSIS-PTS', 'EMP-GLYCOLYSIS-PTS-TCA'}; % 'TCA', 
  script_fsc_obd_pathways_display;

  data_id   = 'd3_sce'; run_id = 'FSC_sce_KcatBrenda'; kcat_usage = 'use'; organism='sce';
  model_names = {'EMP-GLYCOLYSIS', 'EMP-GLYCOLYSIS-TCA'}; %  'TCA', 
  script_fsc_obd_pathways_display;
end

model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';

for it = 1:length(model_names),

  model_name = model_names{it};
  display(sprintf('\n%d: %s',it,model_name));
  filenames = fsc_filenames(model_dir, model_name, data_id, run_id,organism);

  load(filenames.network_file);
  load(filenames.kinetic_data_file);
  load(filenames.fsc_result_file);

  fsc_display(model_name,network,v,fsc_options,c,u,u_cost,up,A_forward,r,kinetic_data);

end
