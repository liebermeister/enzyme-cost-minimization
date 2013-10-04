% --------------------------------------------------------------------------------------
% Run FSC analysis for a single pathway model (EMP-GLYCOLYSIS, from series 'obd_pathways')
% --------------------------------------------------------------------------------------

% HIER WEITER AUFRAEUMEN MIT reaktionsnamen usw

% UEBERSETZUNG O2/ETHANOL-FLUSS siehe jozsefs arbeit, fig 2D
% fuer hefe: lieber glycerol/ethanol als lactate, daten in jozsefs paper fig 4
% ist laktat realistisch?

organism   = 'eco';  
model_name = 'EMP-GLYCOLYSIS-PTS-TCA-OVERFLOW';
run_id     = 'r1_ferm_resp_eco_brenda'; 
data_id    = 'd1_eco'; 

% organism   = 'sce';  
% model_name = 'EMP-GLYCOLYSIS-TCA-OVERFLOW';
% run_id     = 'r1_ferm_resp_sce_brenda'; 
% data_id    = 'd3_sce'; 

kcat_usage = 'use';
model_dir  = '/home/wolfram/projekte/flux_specific_enzyme_cost/models/elad/obd_pathways/';
fsc_scores = {'obdw', 'fsc2', 'fsc2sub', 'fsc3', 'fsc3prod', 'fsc4cmr'}; 

%% 'fsc4smr' nicht realistisch


% --------------------------------------------------------------------------------------
% load model and kinetic data (from prepared files)

filenames = fsc_filenames(model_dir, model_name, data_id, run_id, organism);

load(filenames.network_file);
load(filenames.flux_file);
load(filenames.kinetic_data_file);
load(filenames.metabolic_data_file);

network.external(network_find_metabolites(network,'H2O')) = 1;
K = full(network_analyse(network));

v_resp = 2 * K(:,1);
v_ferm = 2 * K(:,2);
 
conc_min(network_find_metabolites(network,'CoA')) = nan;
conc_max(network_find_metabolites(network,'CoA')) = nan;

% Fix glucose concentration to 12 mM

ll = label_names({'D-Glucose'},network.metabolites);
conc_min(ll) = 12;
conc_max(ll) = 12;

% Fix phosphate concentration 

ll = label_names({'Orthophosphate'},network.metabolites);
conc_min(ll) = 10;
conc_max(ll) = 10;

netgraph_concentrations(network,isfinite(conc_min),[],1);

% figure(1001); netgraph_concentrations(network,network.external,[],1);
% figure(1002); netgraph_fluxes(network_CoHid,v_resp);
% figure(1003); netgraph_fluxes(network_CoHid,v_ferm);


% ---------------------------------------------------------

v_resp = 0.01 * v_resp;
v_ferm = 0.01 * v_ferm;

display(sprintf(' Rescaling the fluxes: median flux: %d mM/s',median(v_resp)))


% ---------------------------------------------------------
% set cofactor concentrations as in data values (where these are available)

display(' Replacing cofactor concentrations');
ll = label_names({'ATP','ADP','NADH', 'Nicotinamide adenine dinucleotide' ,'NAD+', 'Orthophosphate','NADPH','NADP+'},network.metabolites);
ll = ll(find(ll));
ll = ll(find(isfinite(c_data(ll))));
conc_min(ll) = c_data(ll);
conc_max(ll) = c_data(ll);


% ---------------------------------------------------------
% prepare options 

% fsc_default options ...

clear fsc_options
fsc_options.model_name         = model_name; 
fsc_options.run_id             = run_id;
fsc_options.fsc_scores         = fsc_scores;
fsc_options.conc_min_default   = 10^-3;
fsc_options.conc_max_default   = 10^2;
fsc_options.conc_min           = conc_min;
fsc_options.conc_max           = conc_max;
fsc_options.lambda_regularisation = 10^-5;
fsc_options.kinetic_data       = kinetic_data;
fsc_options.kcat_usage         = kcat_usage;
fsc_options.kcat_prior_median  = 200; % similar to median in glycolysis+tca
fsc_options.kcat_prior_log10_std = 0.5;
fsc_options.c_data             = c_data;
fsc_options.u_data             = u_data;  
fsc_options.ind_scored_enzymes = 1:length(network.actions);
fsc_options.enzyme_cost_weights = network.enzyme_size(fsc_options.ind_scored_enzymes);
fsc_options.show_graphics      = 1;  
fsc_options.show_metabolites   = network.metabolites;
fsc_options.network_CoHid      = network_CoHid;
fsc_options.print_graphics     = 1;
fsc_options.psfile_dir         = filenames.psfile_dir;
fsc_options.quantity_info_file = [filenames.resource_dir '/data-kinetic/quantity_info.tsv'];


% run
% in getrennten modellen kommt fuer pyruvat beide male ugf 0.1 mM raus (bei fsc4cmr) .. dadurch bleiben
% im kombinierten modell die konzentrationen konstant!
% ausserdem: bei fsc4-scores kann es immer noch numerische probleme geben (nicht-konvex!)

alpha_list = 0:0.05:1;

clear my_c my_u;
for it = 1:length(alpha_list),
 alpha =  alpha_list(it)
 [my_c{it},my_u{it}] = flux_specific_enzyme_cost(network, [1-alpha] * v_ferm + [alpha] * v_resp, fsc_options);
end

% -------------------------------------------------------
% grafik

%reaction_order   = load_any_table(filenames.reaction_order_file);
%metabolite_order = load_any_table(filenames.metabolite_order_file);

for itt = 1:length(fsc_scores),
 fsc_score =  fsc_scores{itt};
 clear c_list u_list
 for it = 1:length(alpha_list),
   c_list(:,it) = my_c{it}.(fsc_score);
   u_list(:,it) = my_u{it}.(fsc_score);
 end
 
 figure(10+itt); clf; im(log10(u_list),[-9,-1],network.actions,alpha_list); colorbar; title(fsc_score);
 figure(20+itt); clf; im(log10(c_list),[],network.metabolites,alpha_list); colorbar; title(fsc_score);
 figure(30+itt); clf; 
 subplot(2,2,1); plot(alpha_list,log10(u_list'));  title([ fsc_score ' log10 enzyme levels']);
 subplot(2,2,3); plot(alpha_list,sum(u_list),'.-'); title([ fsc_score ' Enzyme cost']);
 subplot(1,2,2); plot(alpha_list,log10(c_list')); title([ fsc_score ' log10 metabolite levels']); axis tight;
end


save('~/dummi', 'alpha_list', 'my_c', 'my_u', 'fsc_options', 'organism', 'model_name', 'run_id', 'data_id', 'kcat_usage', 'model_dir', 'fsc_scores');
