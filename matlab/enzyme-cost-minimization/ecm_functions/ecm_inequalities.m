% MEASURES_FOR_ENZYME_COSTS_INEQUALITIES - Helper function 
%
% [delta_G_by_RT,eq_cons] = measures_for_enzyme_costs_inequalities(x,N_forward,log_Keq_forward)
%

function [delta_G_by_RT,eq_cons] = measures_for_enzyme_costs_inequalities(x,N_forward,log_Keq_forward)

% numerical accuracy for negativity constraints

epsilon = -10^-15;

delta_G_by_RT = N_forward' * x - log_Keq_forward + epsilon;

eq_cons = [];