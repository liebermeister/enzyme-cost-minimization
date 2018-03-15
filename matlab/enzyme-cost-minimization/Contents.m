% Enzyme Cost Minimization
%
% Main functions 
%   ecm_paramater_balancing      - Run parameter balancing with standard settings
%   ecm_enzyme_cost_minimization - Run enzyme cost minimization
%   ecm_simple                   - Wrapper function for parameter balancing or enzyme cost minimization
% 
% Demo scripts
%   demo_ecm_simple
%   demo_ecm_ecoli_noor_2016
%
% For a list of possible cost scores (minimisation of enzyme and metabolites levels, as well as fitting of 
% enzyme and metabolite data), see ecm_get_score.m
%
% Matlab Toolboxes required
%   Metabolic Network Toolbox    - Metabolic networks    (https://github.com/wolframliebermeister/mnt)
%   SBMLtoolbox                  - SBML import / export  (http://sbml.org/Software/SBMLToolbox)
%   SBtab toolbox                - SBtab format          (https://github.com/wolframliebermeister/sbtab-matlab)
%   efmtool                      - Elementary flux modes (http://www.csb.ethz.ch/tools/efmtool)
%
% (C) 2018 Wolfram Liebermeister  <wolfram.liebermeister@gmail.com>