Build Elad's obd_pathway models and run FSC analysis

1. Where are the data?

General data          files in directory ~/projekte/flux_specific_enzyme_cost/resources/    

Models                files in /home/wolfram/projekte/flux_specific_enzyme_cost/models/elad

 original data file:  original_files/[PATHWAY_SET].txt  
                      (example: obd_pathways.txt manually updated, see original_files/README:)
                      + a copy in [PATHWAY_SET]/[PATHWAY_SET].txt
 sbtab files:         in [PATHWAY_SET]/sbtab
 matlab files         in [PATHWAY_SET]/matlab

 in matlab, filenames are managed by fsc_filenames

2. How does the workflow work?

o pathway2sbtab.py              Extract Models -> SBtab files in subdir [PATHWAY_SET]/sbtab'
                                (see ~/projekte/flux_specific_enzyme_cost/models/elad/README)

o fsc_workflow_prepare_all.m    extract information that concerns all models in the set

o fsc_workflow_prepare_single.m prepare data for each single model

o script_fsc_EMP_GLYCOLYSIS.m   run FSC analysis for glycolysis models (test script)

o script_fsc_obd_pathways.m     run FSC analysis for many models


Run everything for some example models: script_fsc_all