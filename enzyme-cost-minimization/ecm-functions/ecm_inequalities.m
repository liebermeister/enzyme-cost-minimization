% ECM_INEQUALITIES - Helper function
%
% [delta_G_by_RT,eq_cons] = ecm_inequalities(x,N_forward,log_Keq_forward,Theta_min,v)
%
% The flux vector v is optional: if it is omitted, all reactions are treated as active
% (this is used when the arguments have already been restricted to the active reactions)

function [delta_G_by_RT,eq_cons] = ecm_inequalities(x,N_forward,log_Keq_forward,Theta_min,v)

if nargin < 5,
  v = ones(size(log_Keq_forward));
end

ind_v_active = find(v~=0);

delta_G_by_RT = - [log_Keq_forward(ind_v_active)-Theta_min(ind_v_active)-[N_forward(:,ind_v_active)' * x]];

eq_cons = [];
