% Tests for the ECM workflow with the E. coli model (as in demo_ecm_simple)

setup_test_path;
filenames = filenames_demo_ecm_simple;
outdir = tempname;
mkdir(outdir);
cleanup = onCleanup(@() rmdir(outdir,'s'));

%% version number
assert(~isempty(regexp(ecm_version, '^\d+\.\d+\.\d+$', 'once')));

%% resource files are found
assert(exist(filenames.model_data_file,'file') == 2);
assert(exist(filenames.validation_data_file,'file') == 2);

%% parameter balancing
options = struct('actions', 'parameter balancing');
[report, errors] = ecm_simple(filenames.model_data_file, filenames.validation_data_file, outdir, options);
assert(strcmp(errors, 'No errors'));
assert(contains(report, 'Parameter balancing finished'));
assert(exist([outdir filesep 'pb_result_Model.tsv'],'file') == 2);

%% enzyme cost minimization
options = struct('actions', 'ecm', 'ecm_scores', {{'emc4cm'}}, 'replace_cofactors', {{}}, 'fluctuations_safety_margin', 0);
[report, errors] = ecm_simple(filenames.model_data_file, filenames.validation_data_file, outdir, options);
assert(strcmp(errors, 'No errors'));
assert(contains(report, 'ECM finished'));
assert(exist([outdir filesep 'ecm_result_ModelState.tsv'],'file') == 2);
assert(exist([outdir filesep 'ecm_result_StateRuns.tsv'],'file') == 2);
