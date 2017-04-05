% ECM_ONE_RUN - Perform one ECM run
% 
% function [my_c, my_u, my_up, my_u_cost, my_A_forward, my_x, my_grad, my_lambda] = ecm_one_run(ecm_score,pp,x_min,x_max,x_init,lambda_regularisation,opt)
%
% my_up is the same as my_up, except for the fact that the levels of unscored enzymes are replaced by nan
%   in the case of multiple-condition ECM, my_up contains the preemptive protein levels (i.e., only one value per reaction)

function [my_c, my_u, my_up, my_u_cost, my_A_forward, my_x, my_grad, my_lambda] = ecm_one_run(ecm_score,pp,x_min,x_max,x_init,lambda_regularisation,opt)

%% optimize log concentration profile

[my_x, my_fval,my_exitflag,my_output,my_lambda,my_grad] = fmincon(@(xx) ecm_get_score(ecm_score,xx,pp) + ecm_regularisation(xx,x_min,x_max,lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) ecm_inequalities(xx,pp.N_forward,pp.log_Keq_forward),opt);

if my_exitflag < 1, my_exitflag
  error('optimisation unsuccessful'); end 

%if my_exitflag ~= 1, my_exitflag
%  warning('optimisation result may be problematic'); end 

my_c         = exp(my_x);
my_A_forward = RT * [pp.log_Keq_forward - pp.N_forward' * log(my_c)]; my_A_forward(pp.ind_not_scored) = nan;

%% compute resulting enzyme profile

[my_u_cost, my_u, my_u_preemptive] = ecm_get_score(ecm_score, my_x, pp);

my_up                    = my_u;
my_up(pp.ind_not_scored) = nan;

if pp.multiple_conditions, 
  my_up  = my_u_preemptive; 
end
