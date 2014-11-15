function [ineq_constraints, eq_constraints] = ecm_below_threshold(fsc_score,xx,pp,u_threshold,x_min,x_max,fsc_options)

% function needed for enzyme cost minimization

a1 = ecm_get_score(fsc_score,xx,pp) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation) - u_threshold;
 
a2 = fsc_inequalities(xx,pp.N_forward,pp.log_Keq_forward);

ineq_constraints = [a1; a2];
eq_constraints   = [];