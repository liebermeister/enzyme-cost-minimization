% Functions for Enzyme Cost Minimization
%
% Main functions 
%   fsc_paramater_balancing      - Run parameter balancing
%   fsc_enzyme_cost_minimization - Run enzyme cost minimization
% 
% Installation: 
%   Edit fsc_setup.m: insert the path to your data directory
% 
% Demos
%   See directory 'demo'
%
% MATLAB Toolboxes required
%   Metabolic Network Toolbox - (https://github.com/wolframliebermeister/mnt)
%   SBMLtoolbox               - SBML import / export  (see http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox             - SBtab format (https://github.com/wolframliebermeister/sbtab-matlab)
%   Tensor toolbox            - (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html)
%   efmtool                   - Elementary flux modes (see http://www.csb.ethz.ch/tools/efmtool)
%
% (C) 2011-2014
% Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>