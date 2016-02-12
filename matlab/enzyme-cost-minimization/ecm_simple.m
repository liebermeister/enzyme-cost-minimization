function [report, errors] = ecm_simple(infile, outdir, options)

% [report, errors] = ecm_simple(infile, outdir, options)

% This function reads the input file [infile] (in SBtab format), performs 
% Parameter Balancing or ECM, and saves the results to an SBtab file in [outdir]
%
% Possible actions: (in options.actions) 
%
% 'ecm_standard':
%    runs an ECM using standard settings, 
%    input file: prepared model with data ("ModelData")
%
% 'parameter_balancing': 
%    runs an ECM using standard settings, 
%    input file: prepared model with data ("ModelData")
% 
% The function can also be called via the python script 'ecm.py'

errors = '';
report = '';

eval(default('options','struct'));
options_default = struct('actions','ecm_standard','make_report',0);
options         = join_struct(options_default, options);

switch options.actions,
  
  case 'ecm_standard',

    display('ecm_simple: Running ECM');

    outfile = [outdir '/ecm_result.tsv'];
    report  = '';
    errors   = '';

    %% Load model and data from input file
    [my_network, v, c_data, u_data, conc_min, conc_max, my_positions,errors] = ecm_load_model_and_data_sbtab(infile,outdir);
    % if length(errors),
    %   report = sprintf('An error occurred while parsing the SBtab file // No result file saved');
    %   errors = sprintf('An error occurred while parsing the SBtab file: %s', error);
    % end
    errors = '';
    
    if isempty(errors),
      if length(my_positions),
        my_network = netgraph_read_positions(my_network, my_positions,[],1,0,my_network.actions);
      end
      
      %% ECM standard options

      my_ecm_options                      = ecm_default_options(my_network);
      my_ecm_options.conc_min             = conc_min; 
      my_ecm_options.conc_max             = conc_max;
      my_ecm_options.conc_min_default     = 10^-10; 
      my_ecm_options.conc_max_default     = 10^10; 
      my_ecm_options.ecm_scores           = {'ecf1', 'mdf', 'ecf2s', 'ecf2sp', 'ecf3s', 'ecf3sp', 'ecf4geom', 'ecf4cmr'}; % {'mdf'};%
      my_ecm_options.c_data               = c_data;
      my_ecm_options.u_data               = u_data;
      my_ecm_options.Keq_upper            = 100000000;
      my_ecm_options.insert_Keq_from_data = 1;
      my_ecm_options.replace_cofactors    = {'ATP','ADP','Orthophosphate','NADH', 'NAD+', 'NADPH','NADP+','Ubiquinone', 'Ubiquinol'};
      my_ecm_options.compute_tolerance    = 0;
      my_ecm_options.cost_tolerance_factor  = 1.01;  
      my_ecm_options.tolerance_from_hessian = 0;
      my_ecm_options = ecm_update_options(my_network, my_ecm_options, my_network.kinetics);
      
      %% Run ECM      

      %try
        [c,u,u_cost,up,A_forward,mca_info,c_min,c_max,u_min,u_max] = ecm_enzyme_cost_minimization(my_network, my_network.kinetics, v, my_ecm_options);
        ecm_save_result_sbtab(outfile, my_network, c, u, A_forward, struct('flag_one_file',1,'r',my_network.kinetics));
        report = sprintf('ECM finished // Results saved to file %s', outfile);
        errors = '';
      %catch err
      %  display('error!!!!!');
      %  report = sprintf('An error occurred during ECM // No result file saved');
      %  errors  = [ 'An error occurred during ECM: ' err.identifier];
      %end
    end
      
  case 'parameter_balancing',
    
      display('ecm_simple: Running parameter balancing');

      outfile              = [outdir '/pb_result.tsv'];
      options.make_report  = 0;
      options.use_kegg_ids = [1];
      report               = '';
      errors                = '';

      %try
        my_sbtab = sbtab_document_load_from_one(infile);
      %catch err
      %  report = sprintf('An error occurred while parsing the SBtab file // No result file saved');
      %    errors  = sprintf('An error occurred while parsing the SBtab file: %s- Probably the input file is incomplete or syntactically wrong', err.identifier);
      %end

      if isempty(errors),
        try
         errors = network_sbtab_check_for_errors(my_sbtab);
        catch err
          report = sprintf('An error occurred while converting SBtab file to model // No result file saved');
          errors  = sprintf('An error occurred while converting SBtab file to model: %s- Probably the input file is incomplete or syntactically wrong', err.identifier);
        end
      end
      
      if isempty(errors),
        try
         my_network = sbtab_to_network(my_sbtab,struct('load_quantity_table',0));
        catch err
          report = sprintf('An error occurred while converting SBtab file to model // No result file saved');
          errors  = sprintf('An error occurred while converting SBtab file to model: %s- Probably the input file is incomplete or syntactically wrong.', err.identifier);
        end
      end

      if isempty(errors),
        %try
          import_quantity_list = {'standard Gibbs energy of reaction', ...
                              'standard chemical potential','Michaelis constant',...
                              'activation constant', 'inhibitory constant',...
                              'equilibrium constant','substrate catalytic rate constant', ...
                              'product catalytic rate constant'};
          use_kegg_ids  = 1;
          organism_long = [];% 'Escherichia coli';
          sbtab_table_save(my_sbtab.tables.RateConstant,struct('filename','/tmp/my_pb_data.tsv'));
          %% store input rate constants in intermediate file '/tmp/my_pb_data.tsv'
          my_kinetic_data = data_integration_load_kinetic_data(import_quantity_list, [], my_network, '/tmp/my_pb_data.tsv', 0, options.use_kegg_ids, 0, 1, 'Organism',organism_long);
          
          my_ecm_options = ecm_default_options(my_network);
          my_ecm_options.insert_Keq_from_data = 0;
          my_ecm_options.GFE_fixed = 0;
          my_ecm_options.show_graphics = 0;
          
          [r, r_orig, kinetic_data] = ecm_parameter_balancing(my_network, my_ecm_options, my_kinetic_data);
          my_network.kinetics = r;
          my_network.kinetics.type = 'cs';
          network_to_sbtab(my_network,struct('filename',outfile,'write_concentrations',0,'save_in_one_file',1,'write_enzyme_concentrations',0));
          report = sprintf('Parameter balancing finished // Results saved to file %s', outfile);
        %catch err
        %  report = sprintf('An error occurred during parameter balancing // No result file saved');
        %  errors  = [ 'An error occurred during parameter balancing: ' err.identifier];
        %end
      end

end

if isempty(errors),

  if options.make_report,
  display('Making graphics files')
  r                             = my_network.kinetics;
  kinetic_data                  = [];
  my_ecm_options.show_graphics  = 1;
  my_ecm_options.print_graphics = 1;
  my_ecm_options.psfile_dir     = [outdir '/ps-files/' ];
  my_ecm_options.model_id       = 'MODEL';  % for filenames
  my_ecm_options.metabolite_order_file = [];
  my_ecm_options.reaction_order_file   = [];
  eval(sprintf('mkdir %s',my_ecm_options.psfile_dir));
else
  my_ecm_options.show_graphics  = 0;    
end

if my_ecm_options.show_graphics,
  ecm_display(my_ecm_options.model_id,my_network,v,options,my_ecm_options,c,u,u_cost,up,A_forward,r,kinetic_data,c_min,c_max,u_min,u_max);
end
end

if isempty(errors),
  errors = 'No errors';
else
  display(errors)
end
