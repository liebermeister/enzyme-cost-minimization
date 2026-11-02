% ECM_INEQUALITIES - Helper function 
%
% [delta_G_by_RT,eq_cons] = ecm_inequalities(x,N_forward,log_Keq_forward,Theta_min,v)
%

function [delta_G_by_RT,eq_cons] = ecm_inequalities(x,N_forward,log_Keq_forward,Theta_min,v)

ind_v_active = find(v~=0);

delta_G_by_RT = - [log_Keq_forward(ind_v_active)-Theta_min(ind_v_active)-[N_forward(:,ind_v_active)' * x]];

eq_cons = [];