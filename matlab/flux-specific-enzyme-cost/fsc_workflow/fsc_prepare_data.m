function [kinetic_data, c_data, u_data, fsc_options, filenames] = fsc_prepare_data(model_directory, data_id)

% [kinetic_data, c_data, u_data, fsc_options, filenames] = fsc_prepare_data(model_directory, data_id)
%
% Prepare data for FSC model


% TEST

if 0,
  model_directory = '~/projekte/flux_specific_cost/fsc/models/ecoli_glycolysis/';
  data_id = 'brenda'; % 'forward', 'none'
  fsc_prepare_data(model_directory, data_id);
end


% -----------------------------------------------
% read options file 'fsc_options.txt' from model directory

eval(default('data_id','[]'));

my_fsc_options = fsc_read_options([ model_directory '/fsc_options.txt']);
filenames      = fsc_filenames(model_directory, model_directory, my_fsc_options.model_name);
network        = [];

load(filenames.network_file);


% -----------------------------------------------
% read options for data integration

fsc_options.verbose               = 1;
fsc_options.flag_use_kinetic_data = 'brenda';
fsc_options.protein_file          = {};

this_dir  = [fileparts(which(mfilename))];
data_dir  = [this_dir '/../../data/'];

fsc_options.quantity_files = {[ data_dir '/formation_energies.tsv'], [ data_dir '/metabolite_levels_ecoli.tsv'], [ data_dir '/kinetic_constants_brenda.tsv']};

if length(data_id),
  fsc_options = fsc_read_options([ model_directory '/fsc_options_data_' data_id '.txt'], fsc_options);
end  


% -----------------------------------------------
% prepare kinetic data 

if fsc_options.verbose,
  display(sprintf('Generating kinetic data table for model %s / data set %s', filenames.model_name, data_id));
end

kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','Michaelis constant','activation constant',  'inhibitory constant','equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}, [], network,  fsc_options.quantity_files, 0, 1);

switch fsc_options.flag_use_kinetic_data,
  
  case 'brenda',  %% use brenda data
  case 'none',    %% do not use any kinetic constants KM, KA, KI, Kcat

    kk = kinetic_data;
    em = nan * kk.KM.median;
    kk.KM.median  = em; 
    kk.KM.mean    = em; 
    kk.KM.std     = em; 
    kk.KM.mean_ln = em; 
    kk.KM.std_ln  = em;
    em = nan * kk.KM.median;
    kk.KA.median  = em; 
    kk.KA.mean = em; 
    kk.KA.std = em; 
    kk.KA.mean_ln = em; 
    kk.KA.std_ln = em;
    em = nan * kk.KI.median;
    kk.KI.median = em; 
    kk.KI.mean = em; 
    kk.KI.std = em; 
    kk.KI.mean_ln = em; 
    kk.KI.std_ln = em;
    em = nan * kk.Kcatf.median;
    kk.Kcatf.median = em; 
    kk.Kcatf.mean = em; 
    kk.Kcatf.std = em; 
    kk.Kcatf.mean_ln= em; 
    kk.Kcatf.std_ln = em;
    em = nan * kk.Kcatr.median;
    kk.Kcatr.median = em; 
    kk.Kcatr.mean = em; 
    kk.Kcatr.std = em; 
    kk.Kcatr.mean_ln= em; 
    kk.Kcatr.std_ln = em;

    kinetic_data = kk;
    
 case 'forward',  %% make forward catalytic constants similar to 20

   ind_p = find(v>=0);
   ind_m = find(v<0);
   emp   = ones(size(ind_p));
   emm   = ones(size(ind_m));
   kk.Kcatf.median(ind_p)  = 20     * emp; 
   kk.Kcatf.mean(ind_p)    = 20     * emp; 
   kk.Kcatf.std(ind_p)     = 0.1    * emp; 
   kk.Kcatf.mean_ln(ind_p) = log(2) * emp; 
   kk.Kcatf.std_ln(ind_p)  = 0.01   * emp;
   kk.Kcatr.median(ind_m)  = 20     * emm; 
   kk.Kcatr.mean(ind_m)    = 20     * emm; 
   kk.Kcatr.std(ind_m)     = 0.1    * emm; 
   kk.Kcatr.mean_ln(ind_m) = log(2) * emm; 
   kk.Kcatr.std_ln(ind_m)  = 0.01   * emm;

   kinetic_data = kk;

end

% save as .mat file
save(filenames.kinetic_data_file,'kinetic_data');

% save as SBtab file
data_integration_save_kinetic_data(kinetic_data, network, [filenames.kinetic_data_file '.tsv'])


% ---------------------------------------------------------------
% metabolic data 

if fsc_options.verbose,
  display(sprintf('Generating metabolic data table for model %s / data set %s', filenames.model_name, data_id));
end

quantity_info = data_integration_load_quantity_info({'forward enzyme mass action term','reverse enzyme mass action term','pH'}, fsc_options.quantity_info_file);

[network, v, kinetic_data_2, quantity_info, model_quantities, basic_quantities, protein_data] = psa3_pb_prepare(filenames, [], quantity_info, fsc_options.quantity_files, fsc_options.protein_file);

save(filenames.fsc_input_file,'network', 'v', 'kinetic_data_2', 'quantity_info', 'model_quantities', 'basic_quantities', 'protein_data', 'fsc_options');

c_data = kinetic_data_2.c.mean;
if isfield(protein_data,'u'), 
  u_data = protein_data.u.mean;
else
  u_data = [];
end
  
save(filenames.metabolic_data_file,'c_data','u_data');

display(sprintf('Writing data files to directory "%s"',model_directory));
