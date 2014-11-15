function f = fsc_root_dir()

f = [fileparts(which(mfilename)) '/../../../'];
