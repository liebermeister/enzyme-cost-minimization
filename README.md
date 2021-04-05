Enzyme Cost Minimization
------------------------

[Enzyme cost minimization](https://www.metabolic-economics.de/enzyme-cost-minimization/) (ECM) is a method for predicting optimal metabolite and enzyme concentrations in metabolic systems.
This repository contains Matlab code for ECM. Elad Noor provides an [implementation in python](https://gitlab.com/equilibrator/equilibrator-pathway).
If you use ECM in your work, please cite our article *Noor et al. (2016)* (reference below). 

## Matlab dependencies

For some of the MATLAB functions, the following MATLAB toolboxes must be installed
- [SBML toolbox](http://sbml.org/Software/SBMLToolbox) (optional - needed only if SBML files are used)
- [DERIVESTsuite](http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation) (needed only if the option ecm_options.compute_hessian is set)
- Clone the following [GitHub](https://github.com/liebermeister) repositories
    - [`matlab-utils`](https://github.com/liebermeister/matlab-utils) - utility functions
    - [`metabolic-network-toolbox`](https://github.com/liebermeister/metabolic-network-toolbox) - metabolic network toolbox
    - [`sbtab-matlab`](https://github.com/liebermeister/sbtab-matlab) - SBtab toolbox
-  Make sure all the directories and subdirectories are included in your Matlab path

The following packages are optional dependencies of the Metabolic Networks Toolbox, but they are not used in Enzyme Cost Minimization
- [Tensor toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
- [efmtool](http://www.csb.ethz.ch/tools/efmtool)

## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](mailto:wolfram.liebermeister@gmail.com) with any questions or comments.


## References
1. E. Noor, A. Flamholz, A. Bar-Even, D. Davidi, R. Milo, W. Liebermeister (2016), [*The Protein Cost of Metabolic Fluxes: Prediction from Enzymatic Rate Laws and Cost Minimization*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005167), PLOS Comp. Biol., [DOI: 10.1371/journal.pcbi.1005167](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5094713/)
2. M.T. Wortel, E. Noor, M. Ferris, F.J. Bruggeman, W. Liebermeister (2018),
[*Metabolic enzyme cost explains variable trade-offs between microbial growth rate and yield*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006010), PLoS Computational Biology 14(2): e1006010