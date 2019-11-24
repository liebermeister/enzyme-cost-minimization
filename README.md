Enzyme Cost Minimization
==========================================

Enzyme cost minimization (ECM) is a method for predicting optimal metabolite and enzyme concentrations in metabolic systems. For a description and
further references, please see the [ECM website](https://www.metabolic-economics.de/enzyme-cost-minimization/).

This repository contains Matlab code only. Python users, please consider Elad Noor's ECM implementation at [gitlab](https://gitlab.com/elad.noor/enzyme-cost-minimization).

If you use enzyme cost minimization in your work, please cite our article

Noor E., Flamholz A., Bar-Even A., Davidi D., Milo R., Liebermeister W. (2016), The protein cost of metabolic fluxes: prediction from enzymatic rate laws and cost minimization, [PLoS Comp. Biol. 12 (10): e1005167](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006010)

## Dependencies: 

For some of the MATLAB functions, the following MATLAB toolboxes must be installed

  o Matlab utility functions    (https://github.com/liebermeister/matlab-utils)

  o Metabolic Network Toolbox (https://github.com/liebermeister/mnt)

  o SBMLtoolbox               (http://sbml.org/Software/SBMLToolbox)

  o SBtab toolbox             (https://github.com/liebermeister/sbtab-matlab)

  o DERIVESTsuite             (http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation)

Please make sure that these matlab packages are installed in your system and that all these directories and subdirectories are included in your matlab path.

Please note that the following packages are required by some functions in the Metabolic Networks Toolbox, but they are not required for Enzyme Cost Minimization

  o Tensor toolbox (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)

  o efmtool        (http://www.csb.ethz.ch/tools/efmtool)


## License
This package is released under the [GNU General Public License](LICENSE).

## Contact
Please contact [Wolfram Liebermeister](wolfram.liebermeister@gmail.com) with any questions or comments.
