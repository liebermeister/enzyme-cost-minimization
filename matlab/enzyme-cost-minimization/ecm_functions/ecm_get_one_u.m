function my_u_it = ecm_get_one_u(it,ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options)

% function my_u_it = ecm_get_one_u(it,ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options)
%
% function needed for enzyme variability calculation after enzyme cost minimization:
% given log metabolite profile xx, evaluate enzyme levels 
% and return the levels of the it'th enzyme

[my_u_cost, my_u] = ecm_get_score(ecm_score,xx,pp);

my_u_it = my_u(it);