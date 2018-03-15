function ecm_options_save(filename,ecm_options)

% ECM_OPTIONS_SAVE - Save options struct in JSON format
%
% ecm_options_save(filename, ecm_options)
  
fileID = fopen(filename,'w');

fprintf(fileID,jsonencode(ecm_options));

