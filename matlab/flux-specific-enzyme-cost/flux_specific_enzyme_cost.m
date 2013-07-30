function [c, u, u_tot, up, A_forward, r, r_orig, fsc_options] = flux_specific_cost(network,v,fsc_options)

% PATHWAY_SPECIFIC_COST - Compute optimal flux-specific costs for given flux distribution
%
% [c, u, u_tot, up, A_forward, r, r_orig] = flux_specific_cost(network, v, fsc_options)
%
% Input 
%   network       metabolic network structure (as in Metabolic Network Toolbox)
%   v             flux mode
%   fsc_options   options struct (for fields and defaul values, see 'fsc_options_default')
%
% Output
%   c                 metabolite concentrations
%   c.data            data vector (if provided)
%   c.fix             fixed concentration vector
%   c.init            initial solution vector 
%   c.[SCORE]         result from optimising the score [SCORE]
%                     
%   u                 enzyme levels (all enzymes)
%   u.data            data vector
%   u.[SCORE]         result from optimising the score [SCORE]
%
%   up                enzyme levels (only scored enzymes as defined by 
%                                    fsc_scores.ind_scored_enzymes)
%   up.[SCORE]        result from optimising the score [SCORE]
%
%   u_tot             sum of scored enzyme levels 
%   u_tot.data        data value
%   u_tot.[SCORE]     result from optimising the score [SCORE]
%
%   A_forward         reaction affinity (in flux direction)
%   A_forward.init    initial solution vector 
%   A_forward.[SCORE] result from optimising the score [SCORE]
%
%   r                 Kinetic constants (from parameter balancing)
%   r_orig            Original kinetic constants (used as input in parameter balancing)
%
%   fsc_options       options struct, possibly updated


% -------------------------------------------------
% fsc_options

[nm,nr] = size(network.N);

fsc_options_default = struct('model_name', 'model','run_id','run','network_CoHid',network,'conc_min_default',10^-6,'conc_max_default',10^2,'conc_min',nan*ones(nm,1),'conc_max',nan*ones(nm,1),'c_data',[],'u_data',[],'kinetic_data',[],'print_graphics',0,'show_graphics',1,'flag_given_kinetics',0,'enzyme_cost_weights',[],'kcat_usage','use','kcat_standard',[], 'lambda_regularisation', 10^-10, 'initial_choice','obd','quantity_info_file',[]);

fsc_options_default.ind_scored_enzymes = 1:length(network.actions);
fsc_options_default.show_metabolites   = network.metabolites;

fsc_options_default.fsc_scores = {'fsc2'};
fsc_options_default.ind_scored_enzymes = 1:nr;
fsc_options_default.show_metabolites   = network.metabolites;

fsc_options = join_struct(fsc_options_default,fsc_options);

if fsc_options.flag_given_kinetics ==0,
  if isempty(fsc_options.kinetic_data),
    fsc_options.kinetic_data = data_integration_load_kinetic_data({'standard chemical potential','Michaelis constant','activation constant',  'inhibitory constant','equilibrium constant','substrate catalytic rate constant', 'product catalytic rate constant'}, [], network, [], 0, 1);
  end
end

opt = optimset('MaxFunEvals',10^20,'TolX',10^-20,'Display','off','Algorithm','active-set');

model_name         = fsc_options.model_name         ; 
network_CoHid      = fsc_options.network_CoHid      ;
ind_scored_enzymes = fsc_options.ind_scored_enzymes ;
conc_min_default   = fsc_options.conc_min_default           ;
conc_max_default   = fsc_options.conc_max_default           ;
conc_min           = fsc_options.conc_min           ;
conc_max           = fsc_options.conc_max           ;
fsc_scores         = fsc_options.fsc_scores         ;
show_metabolites   = fsc_options.show_metabolites   ;
c_data             = fsc_options.c_data             ;
u_data             = fsc_options.u_data             ;  
kinetic_data       = fsc_options.kinetic_data       ;
print_graphics     = fsc_options.print_graphics     ;  
show_graphics      = fsc_options.show_graphics      ;  
if show_graphics == 0, print_graphics = 0; end

ind_conc_fix = find(fsc_options.conc_min == fsc_options.conc_max);
met_fix      = network.metabolites(ind_conc_fix);
conc_fix     = fsc_options.conc_min(ind_conc_fix);

fsc_options.conc_fix = conc_fix;
fsc_options.met_fix  = met_fix;

if length(fsc_options.enzyme_cost_weights),
  enzyme_cost_weights= fsc_options.enzyme_cost_weights / median(fsc_options.enzyme_cost_weights);
else,
  enzyme_cost_weights= ones(length(ind_scored_enzymes),1);
end

ind_met_fix = label_names(met_fix, network.metabolites);

if setdiff(find(network.external),ind_met_fix),
  display('  Some external metabolites in network will not be treated as fixed');
end


% --------------------------------------------------------------------------------------
% network, v, kinetic_data, ind_scored_enzymes, met_fix, conc_fix


[Mplus, Mminus, Wplus, Wminus, nm, nr] = make_structure_matrices(network.N,network.regulation_matrix,find(network.external));

ind_not_scored     = setdiff(1:nr,ind_scored_enzymes);


% --------------------------------------------------------------------------------------
% true concentrations and protein levels

c.data = c_data;
u.data = u_data;
u_tot.data = nansum(u.data);


% --------------------------------------------------------------------------------------
% run parameter balancing for kinetic constants

if fsc_options.flag_given_kinetics,

  switch network.kinetics.type 
    case {'cs','ms','ds'},
      display('Using kinetic constants found in network.kinetics');
      r.Keq     = network.kinetics.Keq;
      r.KM      = network.kinetics.KM;
      [r.Kcatf, r.Kcatr] = modular_KV_Keq_to_kcat(network.N,network.kinetics);
    otherwise error('Kinetics cannot be handled');
  end
  
else
  
  switch fsc_options.kcat_usage
    case 'use',
    
    case 'none',

      %% do not use any kcat data
      kk = kinetic_data;
      em = nan * kk.KM.median;
      kk.KM.median = em; kk.KM.mean = em; kk.KM.std = em; kk.KM.mean_ln = em; kk.KM.std_ln = em;
      em = nan * kk.KM.median;
      kk.KA.median = em; kk.KA.mean = em; kk.KA.std = em; kk.KA.mean_ln = em; kk.KA.std_ln = em;
      em = nan * kk.KI.median;
      kk.KI.median = em; kk.KI.mean = em; kk.KI.std = em; kk.KI.mean_ln = em; kk.KI.std_ln = em;
      em = nan * kk.Kcatf.median;
      kk.Kcatf.median = em; kk.Kcatf.mean = em; kk.Kcatf.std = em; kk.Kcatf.mean_ln= em; kk.Kcatf.std_ln = em;
      em = nan * kk.Kcatr.median;
      kk.Kcatr.median = em; kk.Kcatr.mean = em; kk.Kcatr.std = em; kk.Kcatr.mean_ln= em; kk.Kcatr.std_ln = em;
      kinetic_data = kk;
      
    case 'forward',
      %% invent kcat data such that forward values (along flux) ~= 
      %% match the standard value
      ind_p = find(v>=0);
      ind_m = find(v<0);
      emp   = ones(size(ind_p));
      emm   = ones(size(ind_m));
      
      if isempty(fsc_options.kcat_standard), error('Kcat standard value missing'); end
      kcat_forward_value      = fsc_options.kcat_standard; % unit: 1/s
      kk                      = kinetic_data;
      kk.Kcatf.median(ind_p)  =     kcat_forward_value  * emp; 
      kk.Kcatr.median(ind_m)  =     kcat_forward_value  * emm; 
      kk.Kcatf.mean_ln(ind_p) = log(kcat_forward_value) * emp; 
      kk.Kcatr.mean_ln(ind_m) = log(kcat_forward_value) * emm; 
      kk.Kcatf.std_ln(ind_p)  = 0.1 * emp;
      kk.Kcatr.std_ln(ind_m)  = 0.1 * emm;
      [kk.Kcatf.mean(ind_p),kk.Kcatf.std(ind_p)] = lognormal_log2normal(kk.Kcatf.mean_ln(ind_p),kk.Kcatf.std_ln(ind_p));
      [kk.Kcatr.mean(ind_m),kk.Kcatr.std(ind_m)] = lognormal_log2normal(kk.Kcatr.mean_ln(ind_m),kk.Kcatr.std_ln(ind_m));
      kk.Kcatf.mean(ind_m) = nan;
      kk.Kcatr.mean(ind_p) = nan;
      kk.Kcatf.std(ind_m)  = nan;
      kk.Kcatr.std(ind_p)  = nan;
      kinetic_data         = kk;
  end
    
  [r,r_orig,kinetic_data] = parameter_balancing_kinetic(network, kinetic_data,[],[],struct('kcat_median',fsc_options.kcat_standard,'GFE_fixed',1,'quantity_info_file',fsc_options.quantity_info_file));

  %print_matrix([r_orig.Kcatf, r_orig.Kcatr, r.Kcatf,  r.Kcatr, r.Keq],network.actions);
  network.kinetics = set_kinetics(network,'cs',r);

end


% --------------------------------------------------------------------------------------
% initialise

conc_min(isnan(conc_min)) = conc_min_default;
conc_max(isnan(conc_max)) = conc_max_default;

x_min = log(conc_min);
x_max = log(conc_max);
x_fix = nan * ones(nm,1);
x_fix(ind_met_fix) = log(conc_fix); %  , 
ind_fix = find(isfinite(x_fix)); 

x_min(ind_fix) = x_fix(ind_fix);% - 10^-4;
x_max(ind_fix) = x_fix(ind_fix);% + 10^-4;

ind_pos = find(v>=0);
ind_neg = find(v<0);

M_forward(ind_pos,:)        = Mplus(ind_pos,:);
M_forward(ind_neg,:)        = Mminus(ind_neg,:);
Mprod_forward(ind_pos,:)    = Mminus(ind_pos,:);
Mprod_forward(ind_neg,:)    = Mplus(ind_neg,:);
N_forward(:,ind_pos)        =  network.N(:,ind_pos);
N_forward(:,ind_neg)        = -network.N(:,ind_neg);
kc_forward(ind_pos,1)       = r.Kcatf(ind_pos);
kc_forward(ind_neg,1)       = r.Kcatr(ind_neg);
log_Keq_forward(ind_pos,1)  = log(r.Keq(ind_pos));
log_Keq_forward(ind_neg,1)  = -log(r.Keq(ind_neg));
KM                          = full(r.KM);
KM(network.N'==0)           = 1;
KMf                         = exp(sum(Mplus.*log(KM),2));
KMr                         = exp(sum(Mminus.*log(KM),2));
km_forward(ind_pos,1)       = KMf(ind_pos);
km_forward(ind_neg,1)       = KMr(ind_neg);
kmprod_forward(ind_pos,1)   = KMr(ind_pos);
kmprod_forward(ind_neg,1)   = KMf(ind_neg);


% --------------------------------------------------------------------------------------
% find feasible initial solution

epsilon = 10^-10;

try
  x_init = find_polytope_centre([],[], N_forward', log_Keq_forward - epsilon, x_min, x_max);
catch
  error('*** The concentration constraints seem to be infeasible ***');
  c = [];  u = [];  u_tot = [];  up = [];  A_forward = [];  r = [];  r_orig = [];
  return
end

switch fsc_options.initial_choice,
  case 'interval_center',
  display(' Choosing initial concentrations close to centers of allowed intervals'); 
  x_init = fmincon(@(xx) fsc_regularisation(xx,x_min,x_max,1), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
  case 'obd',
    x_init = fmincon(@(xx) obd(xx,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights) + fsc_regularisation(xx,x_min,x_max,fsc_options.lambda_regularisation), x_init,[],[],[],[],x_min,x_max,@(xx) fsc_inequalities(xx,N_forward,log_Keq_forward),opt);
    display(' Choosing OBD solution as initial concentrations'); 
  case 'polytope_center',
  display(' Choosing initial concentrations at polytope center'); 
end

c_init = exp(x_init);
A_forward_init = RT * [log_Keq_forward  - N_forward' * log(c_init)]; 
A_forward_init(ind_not_scored) = nan;

c.fix          = exp(x_fix);
c.init         = c_init;
A_forward.init = A_forward_init;


% --------------------------------------------------------------------------------------
% run optimisations

% gp  = struct('arrowsize',0.05,'actstyle','none','showsign',0,'actprintnames',0,'flag_edges',1);

for it_method = 1:length(fsc_scores),

  fsc_score = fsc_scores{it_method};

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

  my_c = exp(my_x);

  my_A_forward = RT * [log_Keq_forward  - N_forward' * log(my_c)]; my_A_forward(ind_not_scored) = nan;
  
  %% compute resulting necessary enzyme profile
  
  switch fsc_score,
    
    case {'obd'},
      [my_u_tot, my_u] = obd(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'obdw'},
      [my_u_tot, my_u] = obdw(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'mfsc2sub'},
      [my_u_tot, my_u] = mfsc2sub(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc1'},
      [my_u_tot, my_u] = fsc1(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc2sub'},
      [my_u_tot, my_u] = fsc2sub(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc2'},
      [my_u_tot, my_u] = fsc2(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'mfsc2'},
      [my_u_tot, my_u] = mfsc2(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc3'},
      [my_u_tot, my_u] = fsc3(my_x,v,M_forward,N_forward,log_Keq_forward,kc_forward,km_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc3prod'},
      [my_u_tot, my_u] = fsc3prod(my_x,v,M_forward,Mprod_forward,N_forward,log_Keq_forward,kc_forward,km_forward,kmprod_forward,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4dmr'},
      [my_u_tot, my_u] = fsc4dmr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4smr'},
      [my_u_tot, my_u] = fsc4smr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    case {'fsc4cmr'},
      [my_u_tot, my_u] = fsc4cmr(my_x,v,network,ind_scored_enzymes,enzyme_cost_weights);
    otherwise,
      my_u_tot = nan; my_u= nan * ones(nr,1); my_up = nan * ones(nr,1);
  end

  my_up = my_u; my_up(ind_not_scored) = nan;

  c.(fsc_score)     = my_c;
  u.(fsc_score)     = my_u;
  up.(fsc_score)    = my_up;
  u_tot.(fsc_score) = my_u_tot;
  A_forward.(fsc_score) = my_A_forward;

end
