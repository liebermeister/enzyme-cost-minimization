function gams_save_efm_file(efm,reaction_ids, efm_ids, efm_file_out)

% efm_file_out = '/home/wolfram/projekte/flux_cost_functions/paper_with_meike/matlab/data/TEST_efm_file.csv';
%  gams_save_efm_file(efm,reaction_ids, efm_ids, efm_file_out)

T = [{'efms'},num2cell(reaction_ids'); ...
      num2cell(efm_ids), num2cell(efm)' ];

mytable(T,'comma',efm_file_out);
