function [efm, reaction_ids, efm_ids] = gams_load_efm_file(efm_file)

% efm_file = '/home/wolfram/projekte/flux_cost_functions/paper_with_meike/carlsson_model_gams_files/network3-50/efms.csv';
% [efm,reaction_ids, efm_ids] = gams_load_efm_file(efm_file)

ff = load_any_table(efm_file,',');
reaction_ids = ff(1,2:end)';
efm_ids   = cell_string2num(ff(2:end,1));
efm       = cell_string2num(ff(2:end,2:end))';
