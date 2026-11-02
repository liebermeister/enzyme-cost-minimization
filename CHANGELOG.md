# Changelog

Version numbers have the form MAJOR.MINOR.PATCH. MAJOR changes when existing calls may stop working.

## 1.0.0 (not released yet)

First numbered version.

- Example files with `e_coli` or `ecoli` in their names renamed to `escherichia-coli` (`escherichia_coli` in Matlab function and script names)
- Version function `ecm_version`
- `ecm_BASEDIR` and `ecm_RESOURCEDIR` renamed to `ecm_basedir` and `ecm_resourcedir`; they return the main directory and the `resources` directory, without a final file separator
- Folders `resources/model-files` and `resources/data-files` renamed to `resources/models` and `resources/data`
- Installation notes moved from `INSTALLATION` and the inner `README.md` to the main `README.md`; new `GUIDE.md`
- `CITATION.cff` and `.gitignore` added
- The main functions check the dependencies automatically, once per Matlab session
- New folder `tests` with tests (run with `runtests('tests')`)
- `load_model_and_data_sbtab` and `save_model_and_data_sbtab` moved here from Parameter balancing for Matlab. The validation data file is now optional, and `write_delta_G0` is a real argument of `save_model_and_data_sbtab`
- `ecm_save_result_sbtab` accepts fluxes given as a vector
