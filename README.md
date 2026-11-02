Enzyme Cost Minimization
------------------------

[Enzyme cost minimization](https://www.metabolic-economics.de/enzyme-cost-minimization/) (ECM) is a method for predicting optimal metabolite and enzyme concentrations in metabolic systems.
This repository contains Matlab code for ECM. Elad Noor provides an [implementation in python](https://gitlab.com/equilibrator/equilibrator-pathway).
If you use ECM in your work, please cite our article *Noor et al. (2016)* (reference below). 

## Quick start

After the installation (see below), run one of the demo scripts in Matlab:

```matlab
demo_ecm_simple                      % ECM for an E. coli model, using the wrapper function ecm_simple
demo_ecm_escherichia_coli_noor_2016  % the same model, with the workflow in more detail
```

Results are written to Matlab's `tempdir`. To apply ECM to your own model, adjust these scripts. For a description of the workflow, see the
[ECM website](https://www.metabolic-economics.de/enzyme-cost-minimization/workflow.html); for common tasks, see [GUIDE.md](GUIDE.md).

## Installation

The folder `enzyme-cost-minimization` of this repository, with its subfolders, must be in your Matlab path.
The main functions check once per Matlab session whether the required packages are installed (`ecm_dependencies`).

For some of the MATLAB functions, the following MATLAB toolboxes must be installed
For some of the MATLAB functions, the following MATLAB toolboxes must be installed
- [SBML toolbox](http://sbml.org/Software/SBMLToolbox) (optional - needed only if SBML files are used)
- [DERIVESTsuite](http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation) (needed only if the option ecm_options.compute_hessian is set)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - Utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - Metabolic Network Toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
-  Make sure all the directories and subdirectories are included in your Matlab path

The following packages are optional dependencies of the Metabolic Networks Toolbox, but they are not used in Enzyme Cost Minimization
- [Tensor toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
- [efmtool](http://www.csb.ethz.ch/tools/efmtool)

The code was tested with Matlab R2025b on Linux.

## Contents

| Folder | Contents |
|---|---|
| `enzyme-cost-minimization` | main functions `ecm_simple`, `ecm_enzyme_cost_minimization`, `ecm_parameter_balancing` |
| `enzyme-cost-minimization/ecm-functions` | other ECM functions, with the cost scores in `ecm_scores` |
| `enzyme-cost-minimization/demo` | demo scripts |
| `resources` | example models (`resources/models`) and data (`resources/data`) |
| `tests` | tests |

Type `help enzyme-cost-minimization` for an overview, and `help ecm_default_options` for all options.

## Tests

To check your installation, run in the main directory of this repository:
```matlab
results = runtests('tests')
```

## Version

The version number is given by `ecm_version`; changes are listed in [CHANGELOG.md](CHANGELOG.md).

## License
This package is released under the [GNU General Public License](LICENSE), version 2 or later.

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.

## How to cite
If you use ECM, please cite Noor et al. (2016), reference 1 below. Citation information is also provided in the file `CITATION.cff`.

## References
1. E. Noor, A. Flamholz, A. Bar-Even, D. Davidi, R. Milo, W. Liebermeister (2016), [*The Protein Cost of Metabolic Fluxes: Prediction from Enzymatic Rate Laws and Cost Minimization*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005167), PLOS Comp. Biol., [DOI: 10.1371/journal.pcbi.1005167](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5094713/)
2. M.T. Wortel, E. Noor, M. Ferris, F.J. Bruggeman, W. Liebermeister (2018),
[*Metabolic enzyme cost explains variable trade-offs between microbial growth rate and yield*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006010), PLoS Computational Biology 14(2): e1006010