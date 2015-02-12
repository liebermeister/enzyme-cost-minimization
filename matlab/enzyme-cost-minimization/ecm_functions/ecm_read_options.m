function options = ecm_read_options(options_file, default_options)

eval(default('default_options','struct'));

options = struct('model_name','');
fid     = fopen(options_file);
c       = fgetl(fid);

while c ~= -1,
  eval(c);
  c = fgetl(fid);
end

options = join_struct(default_options, options);
