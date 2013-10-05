function fsc_file_locations = fsc_file_locations()

this_dir  = [fileparts(which(mfilename))];

%data_dir = [this_dir '../data/'];
data_dir = fsc_data_dir;

fsc_file_locations                     = struct;
fsc_file_locations.kinetic_data_infile = {[data_dir '/GFEformation_flamholz_2011_millimolar.tsv'], [data_dir '/kinetic_constants_brenda.tsv'], [data_dir '/kinetic_constants_Ecoli_Scerevisiae.csv']};
fsc_file_locations.quantity_files      = {[data_dir '/concentration_ecoli_Orthophosphate.tsv'], [data_dir '/concentration_ecoli_glucose.tsv']};
%fsc_file_locations.protein_file       = {[filenames.model_dir '/' filenames.model_name '_matched_protein_mRNA.tsv']};
fsc_file_locations.quantity_info_file  = '/home/wolfram/matlab/wolf_packages/mnt_1.1/mnt/data_integration/quantity_info.tsv';
fsc_file_locations.kegg_name_conversion_file = {[data_dir '/name-conversion/kegg_compound_names.tsv']};
fsc_file_locations.flag_use_kinetic_data = 'brenda';
fsc_file_locations.verbose               = 1;
fsc_file_locations.data_dir              = data_dir;
