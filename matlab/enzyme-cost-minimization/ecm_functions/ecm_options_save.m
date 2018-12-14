function ecm_options_save(filename,ecm_options)

% ECM_OPTIONS_SAVE - Save options struct to JSON file
%
% ecm_options_save(filename, ecm_options)

% remove all fields that are structs (because they cannot be exported)
% remove all fields that are strings containing "/" characters, such as filenames (because they cannot be exported)

fn = fieldnames(ecm_options);
for it = 1:length(fn),
  if isstruct(ecm_options.(fn{it})), 
    ecm_options = rmfield(ecm_options,fn{it});
  elseif isstr(ecm_options.(fn{it})), 
    if findstr(ecm_options.(fn{it}),'/'),
      ecm_options = rmfield(ecm_options,fn{it});
    end
  end
end
  
fileID = fopen(filename,'w');

fprintf(fileID,jsonencode(ecm_options));
