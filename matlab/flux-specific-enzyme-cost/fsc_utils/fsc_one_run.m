function [my_c, my_u, my_up, my_u_cost, my_A_forward] = fsc_one_run(fsc_score,v,M_forward,Mprod_forward,kmprod_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,ind_not_scored,enzyme_cost_weights,x_min,x_max,x_init,network,fsc_options,opt)

  %% compute optimal log concentration profile
  
  switch fsc_score,
    case 'obd',
      my_x = fmincon(@(xx) obd(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'obdw',
      my_x = fmincon(@(xx) obdw(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'mfsc2sub',
      my_x = fmincon(@(xx) mfsc2sub(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)  + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc1',
      my_x = fmincon(@(xx) fsc1(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc2sub',
      my_x = fmincon(@(xx) fsc2sub(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)  + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc2',
      my_x = fmincon(@(xx) fsc2(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'mfsc2',
      my_x = fmincon(@(xx) mfsc2(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc3',
      my_x = fmincon(@(xx) fsc3(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights)  + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc3prod',
      my_x = fmincon(@(xx) fsc3prod(xx,v,M_forward,Mprod_forward,N_forward,log_Keq_forward,kc_forward,km_forward,kmprod_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc4dmr',
      my_x = fmincon(@(xx) fsc4dmr(xx,v,network,ind_scored_enzymes,enzyme_cost_weights), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc4smr',
      my_x = fmincon(@(xx) fsc4smr(xx,v,network,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    case 'fsc4cmr',
      my_x = fmincon(@(xx) fsc4cmr(xx,v,network,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
  end

  my_c         = exp(my_x);
  my_A_forward = RT * [log_Keq_forward  - N_forward' * log(my_c)]; my_A_forward(ind_not_scored) = nan;
  
  %% compute resulting necessary enzyme profile
  
  switch fsc_score,
    case {'obd'},
      [my_u_cost, my_u] = obd(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'obdw'},
      [my_u_cost, my_u] = obdw(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'mfsc2sub'},
      [my_u_cost, my_u] = mfsc2sub(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc1'},
      [my_u_cost, my_u] = fsc1(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc2sub'},
      [my_u_cost, my_u] = fsc2sub(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc2'},
      [my_u_cost, my_u] = fsc2(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'mfsc2'},
      [my_u_cost, my_u] = mfsc2(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc3'},
      [my_u_cost, my_u] = fsc3(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc3prod'},
      [my_u_cost, my_u] = fsc3prod(my_x,v,M_forward,Mprod_forward,N_forward,log_Keq_forward,kc_forward,km_forward,kmprod_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4dmr'},
      [my_u_cost, my_u] = fsc4dmr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4smr'},
      [my_u_cost, my_u] = fsc4smr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4cmr'},
      [my_u_cost, my_u] = fsc4cmr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    otherwise,
      my_u_cost = nan; my_u = nan * ones(nr,1);
  end

  my_up                 = my_u; 
  my_up(ind_not_scored) = nan;
