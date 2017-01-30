function [c, u, u_cost, up, A_forward, mca_info, c_min, c_max, u_min, u_max, r, u_capacity, eta_energetic, eta_saturation] = ecm_enzyme_cost_minimization(network,r,v,ecm_options)

% ECM_ENZYME_COST_MINIMIZATION - Compute optimal flux-specific enzyme costs for given flux distribution
%
% [c, u, u_cost, up, A_forward, mca_info, c_min, c_max, u_min, u_max, r, u_capacity, eta_energetic, eta_saturation] = ecm_enzyme_cost_minimization(network, e, v, ecm_options)
%
% Input 
%   network       metabolic network structure (as in Metabolic Network Toolbox)
%   r             Kinetic constants (from parameter balancing)
%   v             flux mode
%   ecm_options   options struct (for fields and default values, see 'ecm_default_options')
%
% Output
%   c                 metabolite concentrations
%   c.data            data vector (if provided)
%   c.fixed           fixed concentration vector
%   c.initial         initial solution vector 
%   c.[SCORE]         1st column: result from optimising the score [SCORE]
%                     other columns: possible sampled solutions
%                     
%   u                 enzyme levels (all enzymes)
%   u.data            data vector
%   u.[SCORE]         1st column: result from optimising the score [SCORE]
%                     other columns: possible sampled solutions
%
%   up                enzyme levels (only scored enzymes as defined by 
%                                    ecm_scores.ind_scored_enzymes)
%   up.[SCORE]        1st column: result from optimising the score [SCORE]
%                     other columns: possible sampled solutions
%
%   u_cost            sum of scored enzyme levels 
%   u_cost.data       data value
%   u_cost.[SCORE]    1st column: result from optimising the score [SCORE]
%                     other columns: possible sampled solutions
%
%   A_forward         reaction affinity (in flux direction)
%   A_forward.initial    initial solution vector 
%   A_forward.[SCORE] 1st column: result from optimising the score [SCORE]
%                     other columns: possible sampled solutions


% --------------------------------------------------------------------------------------
% Enzyme cost minimization
% --------------------------------------------------------------------------------------

opt = optimset('MaxFunEvals',10^15,'MaxIter',10^15,'TolX',10^-10,'Display','off','Algorithm','sqp');

network.kinetics = set_kinetics(network,'cs',r);


% --------------------------------------------------------------------------------------
% initialise 

[Mplus, Mminus, Wplus, Wminus, nm, nr] = make_structure_matrices(network.N,network.regulation_matrix,find(network.external));

ind_scored_enzymes = ecm_options.ind_scored_enzymes ;
ind_not_scored     = setdiff(1:nr,ind_scored_enzymes);

c_data              = ecm_options.c_data             ;
u_data              = ecm_options.u_data             ;  
conc_min_default    = ecm_options.conc_min_default   ;
conc_max_default    = ecm_options.conc_max_default   ;
conc_fix            = ecm_options.conc_fix;
conc_min            = ecm_options.conc_min;
conc_max            = ecm_options.conc_max;
conc_min(isnan(conc_min)) = conc_min_default;
conc_max(isnan(conc_max)) = conc_max_default;
ecm_scores          = ecm_options.ecm_scores         ;
ind_met_fix         = ecm_options.ind_met_fix;
enzyme_cost_weights = ecm_options.enzyme_cost_weights;


% --------------------------------------------------------------------------------------
% initialise c, u, u_cost; insert true concentrations and protein levels

c.data      = c_data;
u.data      = u_data;
u_cost.data = nansum(u.data);


% --------------------------------------------------------------------------------------
% adjust some network- and kinetics-related variables

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


% --------------------------------------------------------------------------
% struct with parameters used in cost score functions

pp.v                   = v                  ; 
pp.network             = network            ;
pp.M_forward           = M_forward          ;
pp.Mprod_forward       = Mprod_forward      ;
pp.N_forward           = N_forward          ;
pp.log_Keq_forward     = log_Keq_forward    ;
pp.kc_forward          = kc_forward         ;
pp.km_forward          = km_forward         ;
pp.kmprod_forward      = kmprod_forward     ;
pp.ind_scored_enzymes  = ind_scored_enzymes ;
pp.enzyme_cost_weights = enzyme_cost_weights;
pp.ind_not_scored      = ind_not_scored     ;


% --------------------------------------------------------------------------------------
% find feasible initial solution (as a test for feasibility) and extreme metabolite profiles

%epsilon = 10^-10;   % flux directions must be possible at least
%epsilon = 1 * 1/RT;  % minimal reaction GFE of 1 kJ/mol is required  
epsilon = 10^-10 * double(v~=0) * 1/RT;  % minimal reaction GFE of 10^-10 kJ/mol is required;

try
  [x_start, x1, x2] = find_polytope_centre([],[], N_forward', log_Keq_forward - epsilon, x_min, x_max, 0*x_min);
catch
  error(sprintf('*** No polytope centre found in search for initial solution (required driving forces > %f RT) - maybe the concentration constraints are too tight? ***\n',epsilon));
  c = [];  u = [];  u_cost = [];  up = [];  A_forward = [];
  return
end

X_extreme = [x1, x2]; % matrix with extreme initial concentration vectors


% --------------------------------------------------------------------------
% choose the starting point

switch ecm_options.initial_choice,
  case 'polytope_center',
    display('  Choosing polytope center as starting point'); 
  case 'interval_center',
    display('  Choosing point close to center of allowed intervals as starting point'); 
    x_start = fmincon(@(xx) ecm_regularisation(xx,x_min,x_max,1), x_start,[],[],[],[],x_min,x_max,@(xx) ecm_inequalities(xx,N_forward,log_Keq_forward),opt);
  case 'mdf',
    x_start = fmincon(@(xx) ecm_mdf(xx,pp) + ecm_regularisation(xx,x_min,x_max,ecm_options.lambda_regularisation), x_start,[],[],[],[],x_min,x_max,@(xx) ecm_inequalities(xx,N_forward,log_Keq_forward),opt);
    display('  Choosing MDF solution as starting point'); 
end


% --------------------------------------------------------------------------------------
% run optimisations starting from x_start

x_init = x_start;

c_init         = exp(x_init);
A_forward_init = RT * [log_Keq_forward  - N_forward' * log(c_init)]; 
A_forward_init(ind_not_scored) = nan;
c.fixed          = exp(x_fix);
c.initial         = c_init;
A_forward.initial = A_forward_init;

display('  Running optimisation using normal starting point');

for it_method = 1:length(ecm_scores),

  ecm_score = ecm_scores{it_method};
  display(sprintf('   %s',ecm_score));
  
  [my_c, my_u, my_up, my_u_cost, my_A_forward, my_x, my_grad, my_lambda] = ecm_one_run(ecm_score,pp,x_min,x_max,x_init,ecm_options,opt);

  c.(ecm_score)         = my_c;
  u.(ecm_score)         = my_u;
  up.(ecm_score)        = my_up;
  u_cost.(ecm_score)    = my_u_cost;
  A_forward.(ecm_score) = my_A_forward;
  mca_info.(ecm_score).x_gradient = my_grad;
  
  if ecm_options.compute_hessian, 
    my_x = log(c.(ecm_score));
    my_fct = @(xx) ecm_get_score(ecm_score,xx,pp) + ecm_regularisation(xx,x_min,x_max,ecm_options.lambda_regularisation);
    mca_info.(ecm_score).x_hessian = hessian(my_fct,my_x);
    mca_info.(ecm_score).rates     = ecm_get_specific_rates(ecm_score,my_x,pp);
    my_fct_2 = @(xx) log(ecm_get_score(ecm_score,xx,pp) + ecm_regularisation(xx,x_min,x_max,ecm_options.lambda_regularisation));
    log_hessian = hessian(my_fct_2,my_x);
    %% fix negative eigenvalues if they occur
    my_eigs = eig(log_hessian);
    if min(my_eigs)<0, log_hessian = log_hessian - 1.001 * min(my_eigs) * eye(nm);end
    mca_info.(ecm_score).x_log_hessian = log_hessian;

  end

  if ecm_options.compute_elasticities,
    %% elasticities between the specific rate (i.e. v/u) and x
    for it_r = 1:nr,
      mca_info.(ecm_score).elasticities_rate_x(it_r,:) = gradest(@(xx) ecm_get_specific_rates(ecm_score,xx,pp,it_r),my_x)';
    end
  end
  
end


% --------------------------------------------------------------------------------------
% run more optimisations starting from extreme points, or with relaxed optimality condition

c_min = [];
c_max = [];
u_min = [];
u_max = [];


for it_method = 1:length(ecm_scores),

  ecm_score = ecm_scores{it_method};
  if ecm_options.multiple_starting_points + ecm_options.compute_tolerance,
    display(sprintf('   %s',ecm_score));
  end

  %% -------------------------------------------------------------------------
  %% run more optimisations starting from extreme points

  if ecm_options.multiple_starting_points,
  
    display(sprintf('     Running optimisation from %d extreme starting points',size(X_extreme,2)));
    
    for itt = 1:size(X_extreme,2),

      %% start from points near the extreme points
      x_init = 0.1 * x_start + 0.9 * X_extreme(:,itt); 
      fprintf('%d ',itt)

      %% same as above
      c_init                  = exp(x_init);
      A_forward_init          = RT * [log_Keq_forward  - N_forward' * log(c_init)]; 
      A_forward_init(ind_not_scored) = nan;
      c.fixed                   = exp(x_fix);
      c.initial(:,1+itt)         = c_init;
      A_forward.initial(:,1+itt) = A_forward_init;
      [my_c, my_u, my_up, my_u_cost, my_A_forward] = ecm_one_run(ecm_score,pp,x_min,x_max,x_init,ecm_options,opt);
      c.(ecm_score)(:,1+itt)         = my_c;
      u.(ecm_score)(:,1+itt)         = my_u;
      up.(ecm_score)(:,1+itt)        = my_up;
      u_cost.(ecm_score)(:,1+itt)    = my_u_cost;
      A_forward.(ecm_score)(:,1+itt) = my_A_forward;
      
    end
    
    %% copy the best result to the first position
    [dum,ind_opt] = min(u_cost.(ecm_score));
    
    c.initial             = [ c.initial(:,ind_opt)                c.initial               ];  
    c.(ecm_score)         = [ c.(ecm_score)(:,ind_opt)         c.(ecm_score)        ];  
    u.(ecm_score)         = [ u.(ecm_score)(:,ind_opt)         u.(ecm_score)        ];  
    up.(ecm_score)        = [ up.(ecm_score)(:,ind_opt)        up.(ecm_score)       ];  
    u_cost.(ecm_score)    = [ u_cost.(ecm_score)(ind_opt)      u_cost.(ecm_score)   ];  
    A_forward.(ecm_score) = [ A_forward.(ecm_score)(:,ind_opt) A_forward.(ecm_score)];
    A_forward.initial     = [ A_forward.initial(:,ind_opt) A_forward.initial];

  end

  %% -------------------------------------------------------------------------
  %% compute tolerance in log concentration profile
  %% (at some maximum allowed percentage increase in cost)

  if ecm_options.compute_tolerance,

    display('    Computing the effects of relaxed optimality assumptions');
  
    x           = log(c.(ecm_score)(:,1));
    u_threshold = ecm_options.cost_tolerance_factor * u_cost.(ecm_score)(1);
  
    my_x_min = x;
    my_x_max = x;
  
    for it = 1:length(x),
      if x(it) ~= x_min(it),
        [my_x, my_x_min(it)] = fmincon(@(xx) xx(it), x,[],[],[],[],x_min,x_max,@(xx) ecm_below_threshold(ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), opt);
      end
      if x(it) ~= x_max(it),
        [my_x, my_x_max(it)] = fmincon(@(xx) -xx(it), x,[],[],[],[],x_min,x_max,@(xx) ecm_below_threshold(ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), opt);
        my_x_max(it) = -my_x_max(it);
      end
    end
  
    c_min.(ecm_score) = exp(my_x_min);
    c_max.(ecm_score) = exp(my_x_max);
  
    %% compute possible variation in enzyme profile u
    %% (at some maximum allowed percentage increase in cost)
    
    my_u     = u.(ecm_score)(:,1);
    my_u_min = my_u;
    my_u_max = my_u;

    if sum(isfinite(my_u)),
      
    for it = 1:length(u.data),
      [my_x, my_u_min(it),exitflag] = fmincon(@(xx) ecm_get_one_u(it,ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), x,[],[],[],[],x_min,x_max,@(xx) ecm_below_threshold(ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), opt);
  
      [my_x, my_u_max(it),exitflag] = fmincon(@(xx) -ecm_get_one_u(it,ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), x,[],[],[],[],x_min,x_max,@(xx) ecm_below_threshold(ecm_score,xx,pp,u_threshold,x_min,x_max,ecm_options), opt);
  
      if ~exitflag, warning(sprintf('Enzyme variability in %s: No valid solution found',network.actions{it})); 
        my_u_max(it) = nan;
      end 
  
      my_u_max(it) = - my_u_max(it);
    end
  
    u_min.(ecm_score) = my_u_min;
    u_max.(ecm_score) = my_u_max;

    end
  end
  
end


% ---------------------------------------------------------------------
% compute efficiency values

kcat_forward      = r.Kcatf;
kcat_forward(v<0) = r.Kcatr(v<0);
u_capacity        = abs(v)./kcat_forward;
eta_energetic     = struct;
eta_saturation    = struct;

for it_method = 1:length(ecm_options.ecm_scores),
  this_ecm_score   = ecm_options.ecm_scores{it_method};
  my_u             = u.(this_ecm_score);
  switch this_ecm_score,
    case {'mdf','emc1'},
      my_eta_energetic = nan * ones(size(v));
      my_eta_energetic(ecm_options.ind_scored_enzymes) = ones(size(ecm_options.ind_scored_enzymes));
    otherwise,
      my_A_forward     = A_forward.(this_ecm_score);
      my_eta_energetic = 1-exp(-my_A_forward/RT);
      my_eta_energetic(my_eta_energetic<0) = 0;
    end
  switch this_ecm_score,
    case {'mdf','emc1','emc2s','emc2sp'},
      my_eta_saturation = nan * ones(size(v));
      my_eta_saturation(ecm_options.ind_scored_enzymes) = ones(size(ecm_options.ind_scored_enzymes));
    otherwise,
      my_eta_saturation = [abs(v) ./ kcat_forward ./ my_eta_energetic] ./ my_u;
  end
  eta_energetic.(this_ecm_score)  = my_eta_energetic ;
  eta_saturation.(this_ecm_score) = my_eta_saturation;
end
