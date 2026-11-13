# Note: standard concentration of equilibrium constants (corrected 2026-11-14)

Files:

* `escherichia-coli_wortel_2018_data.tsv`

The equilibrium constants in these files (E. coli central metabolism, Carlsson / Wortel et al. 2018
model; values from eQuilibrator) refer to the standard concentration **1 M**. The files did not declare
this, and tools read the values as 1 mM. For reactions with different numbers of substrates and
products, this changes Keq by a factor of 1000. In the Wortel model, this concerns R5r (aldolase),
R61r (KDPG aldolase), R54ra and R26r (R26r and R7ra omit phosphate in the model, while the
eQuilibrator values include it).

Correction: the tables now declare `StandardConcentration='1M'` (SBtab table attribute); the values are
unchanged. Results computed earlier from these files with tools that assumed 1 mM are affected.
