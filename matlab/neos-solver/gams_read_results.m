function [metabolites, conc, reactions, enzyme_levels, etot] = gams_read_results(results_dir,model_number, parameter_number, solution_number,nm,nr)

% Read results from gams optimization

C = load_any_table(sprintf('%s/met-m%d-p%d-e%d.csv',results_dir, model_number, parameter_number, solution_number),',');
E = load_any_table(sprintf('%s/enz-m%d-p%d-e%d.csv',results_dir, model_number, parameter_number, solution_number),',');

metabolites = {};
conc        = [];

for it = 1:size(C,1),
  metab       = C{it,1};
  metabolites = [metabolites; {metab(2:end-1)}];
  conc        = [conc; eval(C{it,2})];
end

reactions      = {};
enzyme_levels  = [];

for it = 1:size(E,1),
  rea           = E{it,1};
  reactions     = [reactions; {rea(2:end-1)}];
  enzyme_levels = [enzyme_levels; eval(E{it,2})];
end

etot = sum(enzyme_levels);
