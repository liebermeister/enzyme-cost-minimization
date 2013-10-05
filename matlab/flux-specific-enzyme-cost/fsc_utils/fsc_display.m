function fsc_display(model_name,network,v,fsc_options,c,u,u_tot,up,A_forward,r,kinetic_data)

% PATHWAY_SPECIFIC_COST_DISPLAY - Display results of Pathway Specific Cost Modelling
%
% fsc_display(network,v,fsc_options,c,u,u_tot,up,A_forward,r)

[nm,nr] = size(network.N);

value_colors = rb_colors;

cm            = jet(length(fsc_options.fsc_scores));
ca

% -----------------------------------

reaction_names = network.actions;

if isfield(network,'genes'),
  gene_names     = network.genes;
else
  gene_names     = network.actions;
end
  
if ~isfield(fsc_options,'reaction_order'),
  fsc_options.reaction_order = reaction_names;
else,
  ll = label_names(fsc_options.reaction_order,reaction_names);
  fsc_options.reaction_order = fsc_options.reaction_order(find(ll));
  fsc_options.reaction_order = [column(fsc_options.reaction_order); column(setdiff(reaction_names,fsc_options.reaction_order))];
end

ind_genes        = label_names(reaction_names,fsc_options.reaction_order);
ind_genes_rev    = label_names(fsc_options.reaction_order,reaction_names);
ind_genes_scored = fsc_options.ind_scored_enzymes(ind_genes_rev);

% -----------------------------------

if ~isfield(fsc_options,'metabolite_order'),
  fsc_options.metabolite_order = network.metabolites;
else,
  fsc_options.metabolite_order = [column(fsc_options.metabolite_order); column(setdiff(network.metabolites,fsc_options.metabolite_order))];
end

ind_show_met = [];
for it = 1:length(fsc_options.metabolite_order),
  my_metabolite = intersect(fsc_options.metabolite_order{it},fsc_options.show_metabolites);
  if length(my_metabolite),
    ind_show_met = [ ind_show_met, label_names(my_metabolite, network.metabolites)];
  end
end

% -----------------------------------

gp = struct('arrowsize',0.05,'arrowstyle','fluxes','actstyle','none','arrowvalues',v);
figure(27); netgraph_concentrations(fsc_options.network_CoHid,[],v,1,gp);

figure(4); clf; 
gp = struct('arrowsize',0.05,'actstyle','none','colorbar',1,'showsign',0,'actprintnames',0,'metprintnames',1,'flag_edges',1,'colormap', value_colors);
gp.actvaluesmin = log10(10);
gp.actvaluesmax = log10(10000);
gp.actstyle     = 'fixed';
gp.arrowstyle   = 'none';
netgraph_concentrations(fsc_options.network_CoHid,[],log10(r.Kcatf),1,gp); 
title('log10 Kcat+, balanced (1/s)');

figure(14); clf; 
gp = struct('arrowsize',0.05,'actstyle','none','colorbar',1,'showsign',0,'actprintnames',0,'metprintnames',1,'flag_edges',1,'colormap', value_colors);
gp.actvaluesmin = log10(10);
gp.actvaluesmax = log10(10000);
gp.actstyle     = 'fixed';
gp.arrowstyle   = 'none';
netgraph_concentrations(fsc_options.network_CoHid,[],log10(kinetic_data.Kcatf.median),1,gp); 
title('log10 Kcat+, original (1/s)');


figure(11); clf; 
gp.actstyle = 'none';
gp.metvaluesmin = log10(fsc_options.conc_min_default);
gp.metvaluesmax = log10(100);% (fsc_options.conc_max_default);
netgraph_concentrations(network,log10(c.init(:,1)),[],1,gp); title('Initial state: log10 Concentrations (mM)');

figure(12); clf; 
gp.arrowstyle='none'; gp.actstyle = 'fixed';
gp.actvaluesmin = log10(1);
gp.actvaluesmax = log10(100);
gp.metvaluesmin = [];
gp.metvaluesmax = [];
netgraph_concentrations(fsc_options.network_CoHid,[],log10(A_forward.init),1,gp); title('Inital state: log10 Reaction affinities (kJ/mol)');

figure(8); clf
gp = struct('arrowstyle','none','colorbar',1,'showsign',0,'actprintnames',0,'flag_edges',1,'colormap', value_colors);
gp.metvaluesmin = log10(fsc_options.conc_min_default);
gp.metvaluesmax = log10(100); % fsc_options.conc_max_default);
netgraph_concentrations(network,log10(c.data(:,1)),[],1,gp); title('Experimental data (1st file): log10 Concentrations (mM)');

figure(28); clf
gp = struct('arrowstyle','none','colorbar',1,'showsign',0,'actprintnames',0,'flag_edges',1,'colormap', value_colors);
gp.metvaluesmin = log10(fsc_options.conc_min_default);
gp.metvaluesmax = log10(100); % fsc_options.conc_max_default);
netgraph_concentrations(network,log10(c.fix),[],1,gp); title('Fixed values: log10 Concentrations (mM)');

display('Fixed concentrations')
pm(fsc_options.conc_fix,fsc_options.met_fix)

figure(9); clf
gp = struct('arrowstyle','none','colorbar',1,'showsign',0,'actprintnames',1,'flag_edges',1,'colormap', value_colors);
gp.actnames     = gene_names;
gp.actvaluesmin = log10(nanmin(u.data));
gp.actvaluesmax = log10(nanmax(u.data));
gp.metvaluesmin = [];
gp.metvaluesmax = [];
netgraph_concentrations(fsc_options.network_CoHid,[],log10(u.data),1,gp); title('Experimental data: log10 enzyme levels (a.u.)');


% --------------------------------------
% concentration predictions

% everything in one picture

figure(18); clf; hold on;
set(gcf,'Position',[100 300 1000 500])
xtick = 0:length(ind_show_met)-1;
h     = [];
my_fsc_scores = [{'init'},fsc_options.fsc_scores];
for it = 1:length(my_fsc_scores),
  fsc_score = my_fsc_scores{it};
  h(it)     = plot(xtick,c.(fsc_score)(ind_show_met,1) * exp(0.03*it),'-','MarkerSize',20); %  * exp(0.001*it)
end
line_colors(h,'sunrise_colors');

plot(xtick,c.fix(ind_show_met),'r.','MarkerSize',30);
if length(c.data),
  n_c_data_sets = size(c.data,2);
  plot(repmat(column(xtick),n_c_data_sets,1),reshape(c.data(ind_show_met,:),n_c_data_sets * length(ind_show_met),1),'k.','MarkerSize',20); 
  methods = [my_fsc_scores,{'fixed','data'}];
else
  methods = [my_fsc_scores,{'fixed'}];
end

metnames_show = strrep(network.metabolites(ind_show_met),'D-','');
text(xtick + 0.2,c.data(ind_show_met),metnames_show,'FontSize',10);
legend(methods,'Location','NorthEastOutside','FontSize',10);
axis tight; set(gca,'XTick',xtick,'XTickLabel',metnames_show,'YScale','Log');
my_xticklabel([],fsc_options.conc_min_default,[],10);
title('Comparison: log10 Concentration profiles (mM)');

% only ecf3:

figure(38); clf; hold on;
set(gcf,'Position',[100 300 1000 500])
xtick = 0:length(ind_show_met)-1;
fsc_score = 'fsc3prod';
h    = plot(xtick,c.(fsc_score)(ind_show_met,1),'o','MarkerSize',10);
plot(xtick,c.fix(ind_show_met),'r.','MarkerSize',30);
if length(c.data),
  n_c_data_sets = size(c.data,2);
  plot(repmat(column(xtick),n_c_data_sets,1),reshape(c.data(ind_show_met,:),n_c_data_sets * length(ind_show_met),1),'k.','MarkerSize',20); 
  methods = {fsc_score,'fixed','data'};
else
  methods = {fsc_score,'fixed'};
end

metnames_show = strrep(network.metabolites(ind_show_met),'D-','');
text(xtick + 0.2,c.(fsc_score)(ind_show_met),metnames_show,'FontSize',10);
legend(methods,'Location','NorthEastOutside','FontSize',10);
ylabel('Concentration [mM]');
axis tight; set(gca,'YScale','Log','XTick',[]); 
%'XTick',xtick,'XTickLabel',metnames_show,
a = axis; hold on
for it = 1:length(xtick),
  plot([xtick(it),xtick(it)],[a(3),a(4)],'-');
end
%my_xticklabel([],fsc_options.conc_min_default,[],10);
%title('Comparison: log10 Concentration profiles (mM)');


% --------------------------------------
% enzyme predictions scatter plot

% fsc_score = fsc_options.fsc_scores{end};
fsc_score = 'fsc3prod';

figure(16); clf; hold on; set(gcf,'Position',[150 350 600 600])
ind_finite = find(isfinite(u.data));
my_colors = sunrise_colors(length(u.data));
my_colors = my_colors(ind_finite,:);
for itt = 1:length(ind_finite)
  plot(u.data(ind_finite(itt)),u.(fsc_score)(ind_finite(itt),1),'.','Markersize',30,'Color',my_colors(itt,:)); 
end
text(1.05*u.data(ind_finite),1.05*u.(fsc_score)(ind_finite,1),network.genes(fsc_options.ind_scored_enzymes(ind_finite)),'FontSize',8);
axis tight; axis equal; axis square; 
xlabel('Enzyme level (a.u.)','Fontsize',14);
ylabel('Predicted enzyme level (mM)','Fontsize',14);
set(gca,'XScale','Log','YScale','Log','Fontsize',14);
[cc, pvalue_one_tailed, pvalue_two_tailed] = my_corrcoef(log(u.data(ind_finite)),log(u.(fsc_score)(ind_finite,1)));
title(sprintf('%s R^2: %2.3f // p-value: %2.3f',fsc_score,cc^2,pvalue_one_tailed));


% --------------------------------------
% concentration predictions scatter plot

%fsc_score = fsc_options.fsc_scores{end};
fsc_score = 'fsc3prod';

figure(17); clf; hold on; set(gcf,'Position',[150 350 600 600])

my_c_min = min(nanmin(c.data(:)),nanmin(c.(fsc_score)(:)));
my_c_max = max(nanmax(c.data(:)),nanmax(c.(fsc_score)(:)));

%plot([fsc_options.conc_min_default,fsc_options.conc_max_default],[fsc_options.conc_min_default,fsc_options.conc_max_default],'-');
plot([my_c_min,my_c_max],[my_c_min,my_c_max],'-');

ind_non_fixed = find(label_names(network.metabolites, fsc_options.met_fix)==0);

for ittt = 1:length(ind_non_fixed),
  it = ind_non_fixed(ittt);
  plot(nanmean(c.data(it,:)),c.(fsc_score)(it,1),'k.','Markersize',20);
  plot([nanmin(c.data(it,:)),nanmax(c.data(it,:))],c.(fsc_score)(it,1) * [1 1],'k-');
  text(1.1*nanmean(c.data(it,:)),1.1*c.(fsc_score)(it,1),network.metabolites{it},'FontSize',6);  
end

axis tight; axis equal; axis square; 
xlabel('Median experimental concentration (mM)','Fontsize',18);
ylabel('Predicted concentration (mM)','Fontsize',18);
set(gca,'XScale','Log','YScale','Log','Fontsize',18);

my_nanmean_c_data   = nanmean(c.data,2);
ind_finite          = find(isfinite(my_nanmean_c_data));
[cc, pvalue_one_tailed] = my_corrcoef(log(my_nanmean_c_data(ind_finite)),log(c.(fsc_score)(ind_finite,1)));
title(sprintf('%s R^2: %2.3f // p-value: %2.3f',fsc_score,cc^2,pvalue_one_tailed));


% --------------------------------------
% reaction affinity predictions

figure(19); clf; hold on;
set(gcf,'Position',[200 400 1000 500])
xtick = 0:length(network.actions);
h = [];

my_fsc_scores = [{'init'},fsc_options.fsc_scores];

for it = 1:length(my_fsc_scores),
  fsc_score = my_fsc_scores{it};
  h(it)     = plot(xtick,[0; -cumsum(v .* A_forward.(fsc_score)(:,1))] + (0.002*it),'-'); %  
end
line_colors(h,'sunrise_colors');

legend(my_fsc_scores,'Location','NorthEast','FontSize',16);
axis tight; set(gca,'XTick',xtick,'XTickLabel',[{''};strrep(network.actions,'_', ' ')]);
my_xticklabel([],[],[],10);
title('Cumulative entropy production (v * A)');

% --------------------------------------
% c over KM

figure(20); subplot('Position',[0.15 0.05 0.8 0.7]);
c_over_KM = repmat(c.init(:,1)',length(network.actions),1)./r.KM;
c_over_KM(r.KM==0) = nan;
im(log10(c_over_KM(:,ind_show_met)),[-3,3],network.actions,network.metabolites(ind_show_met)); colormap([0.95 0.95 0.95; value_colors]); colorbar
my_xticklabel
%title('log10 c (initial solution) / KM');

figure(21); subplot('Position',[0.15 0.05 0.8 0.7]);
c_over_KM = repmat(c.fsc3prod(:,1)',length(network.actions),1)./r.KM;
c_over_KM(r.KM==0) = nan;
im(log10(c_over_KM(:,ind_show_met)),[-3,3],network.actions,network.metabolites(ind_show_met)); colormap([0.95 0.95 0.95; value_colors]); colorbar
my_xticklabel
%title('log10 c (initial solution) / KM');


% --------------------------------------
% enzyme level predictions

% barplot for all methods

figure(33); clf; 
if length(u.data),
  ug = u.data(fsc_options.ind_scored_enzymes)/nansum(u.data(fsc_options.ind_scored_enzymes));
  methods = [{'data'},fsc_options.fsc_scores];
else
  ug=[];
  methods = fsc_options.fsc_scores;
end

for it = 1:length(fsc_options.fsc_scores),
  fsc_score = fsc_options.fsc_scores{it};
  ug_meth   = up.(fsc_score)(fsc_options.ind_scored_enzymes,1);
  ug_meth   = ug_meth/sum(ug_meth);
  ug        = [ug, ug_meth];
end
h = bar(flipud(ug(ind_genes_scored,:))','stacked');
colormap(flipud(sunrise_colors)); 
axis([0,4+length(fsc_options.fsc_scores), 0, 1.02]); 
legend(fliplr(h),gene_names(ind_genes_scored),'Fontsize',8,'Location','EastOutside')
set(gca,'XTickLabel',methods); ylabel('Relative protein amount');
set(gca,'FontSize',14);

% color legend

figure(34); clf; 
mmm = nan * ones(nr,1); 
mmm(fsc_options.ind_scored_enzymes) = ind_genes(fsc_options.ind_scored_enzymes);
gpp = struct('actprintnames',1,'arrowstyle','none','showsign',0,'colormap',sunrise_colors,'text_offset',[0.03,0]);
gpp.actnames = gene_names;
netgraph_concentrations(fsc_options.network_CoHid,[],column(mmm),1,gpp); 
title('Enzymes colour legend');

% barplot for fsc3prod only

fsc_score = 'fsc3prod';
ug   = up.(fsc_score)(fsc_options.ind_scored_enzymes,1);
ud   = u.data(fsc_options.ind_scored_enzymes,1);

figure(43); clf; 
cmap = sunrise_colors(length(ind_genes_scored));
subplot(2,1,1);
for it = 1:length(ind_genes_scored);
 h= bar(it,ug(ind_genes_scored(it))); set(h,'FaceColor',cmap(it,:)); hold on
end
subplot(2,1,2);
for it = 1:length(ind_genes_scored);
 h= bar(it,ud(ind_genes_scored(it))); set(h,'FaceColor',cmap(it,:)); hold on
end
set(gca,'XTick',1:length(ind_genes_scored),'XTickLabel',gene_names(ind_genes_scored)); 
ylabel('Predicted enzyme level [mM]');
set(gca,'FontSize',10);
title(sprintf('Total predicted enzyme level (%s): %f mM',fsc_score, sum(ug(ind_genes_scored))));

% --------------------------------------------------------
% plot flux (data) versus enzyme level (data)

figure(44); clf; hold on
ind_finite   = find(isfinite(u.data .* v));
my_colors = sunrise_colors(length(u.data));
my_colors = my_colors(ind_finite,:);
for it = 1:length(ind_finite),
  plot(u.data(ind_finite(it)),v(ind_finite(it)),'.','Color',my_colors(it,:),'MarkerSize',30);
end
text(u.data(ind_finite) * 1.03,v(ind_finite),gene_names(ind_finite));
set(gca,'XScale','Log','YScale','Log');
xlabel('Enzyme level (data)','Fontsize',14); ylabel('Flux (data)','Fontsize',14);
axis equal; axis square; axis tight
[cc, pvalue_one_tailed, pvalue_two_tailed] = my_corrcoef(log(u.data(ind_finite)),log(v(ind_finite)));
title(sprintf('R^2: %2.3f // p-value: %2.3f',cc^2,pvalue_one_tailed));


% plot flux (data) versus enzyme level * kcat (both data)

figure(45);clf; hold on
ind_finite = find(isfinite(u.data .* kinetic_data.Kcatf.median .* v));
my_colors = sunrise_colors(length(u.data));
my_colors = my_colors(ind_finite,:);
for it = 1:length(ind_finite),
  plot(u.data(ind_finite(it)),v(ind_finite(it))/kinetic_data.Kcatf.median(ind_finite(it)),'.','Color',my_colors(it,:),'MarkerSize',30);
end
text(u.data(ind_finite) * 1.03,v(ind_finite)  ./ kinetic_data.Kcatf.median(ind_finite),gene_names(ind_finite));
set(gca,'XScale','Log','YScale','Log');
xlabel('Enzyme level','Fontsize',14); ylabel('Flux (data) / kcat (data)','Fontsize',14);
axis equal; axis square; axis tight
[cc, pvalue_one_tailed, pvalue_two_tailed] = my_corrcoef(log(u.data(ind_finite)),log(v(ind_finite)./kinetic_data.Kcatf.median(ind_finite)));
title(sprintf('R^2: %2.3f // p-value: %2.3f',cc^2,pvalue_one_tailed));


% --------------------------------------------------------
% save graphics

if fsc_options.print_graphics,
  if isfield(fsc_options,'psfile_dir'),
    cd(fsc_options.psfile_dir);
    print([ model_name '_' fsc_options.run_id '_fluxes.eps'],'-f27','-depsc');
    print([ model_name '_' fsc_options.run_id '_data_concentrations.eps'],'-f8','-depsc');
    print([ model_name '_' fsc_options.run_id '_data_enzymes.eps'],'-f9','-depsc');
    print([ model_name '_' fsc_options.run_id '_fix_concentrations.eps'],'-f28','-depsc');
    print([ model_name '_' fsc_options.run_id '_init_concentrations.eps'],'-f11','-depsc');
    print([ model_name '_' fsc_options.run_id '_init_affinities.eps'],'-f12','-depsc');
    print([ model_name '_' fsc_options.run_id '_kcat_balanced.eps'],'-f4','-depsc');
    print([ model_name '_' fsc_options.run_id '_kcat_original.eps'],'-f14','-depsc');
    print([ model_name '_' fsc_options.run_id '_barplot.eps'],'-f33','-depsc');
    print([ model_name '_' fsc_options.run_id '_barplot_legend.eps'],'-f34','-depsc');
    print([ model_name '_' fsc_options.run_id '_FSC_enzyme_scatter.eps'],'-f16','-depsc');
    print([ model_name '_' fsc_options.run_id '_FSC_conc_scatter.eps'],'-f17','-depsc');
    print([ model_name '_' fsc_options.run_id '_all_conc_profiles.eps'],'-f18','-depsc');
    print([ model_name '_' fsc_options.run_id '_affinity_profiles.eps'],'-f19','-depsc');
    print([ model_name '_' fsc_options.run_id '_c_by_KM_initial.eps'],'-f20','-depsc');
    print([ model_name '_' fsc_options.run_id '_c_by_KM_fsc3prod.eps'],'-f21','-depsc');
    print([ model_name '_' fsc_options.run_id '_fsc3prod_conc_prediction.eps'],'-f38','-depsc');
    print([ model_name '_' fsc_options.run_id '_fsc3prod_enzyme_prediction.eps'],'-f43','-depsc');
    print([ model_name '_' fsc_options.run_id '_scatter_enzyme_vs_flux.eps'],'-f44','-depsc');
    print([ model_name '_' fsc_options.run_id '_scatter_enzyme_kcat_vs_flux.eps'],'-f45','-depsc');
  end
end

for it_method = 1:length(fsc_options.fsc_scores),
  fsc_score = fsc_options.fsc_scores{it_method};
  
  my_c = c.(fsc_score)(:,1);
  my_A_forward = A_forward.(fsc_score)(:,1);
  my_up = up.(fsc_score)(:,1);

  figure(1); clf; 
  gp.actstyle = 'none';
  gp.metvaluesmin = log10(fsc_options.conc_min_default);
  gp.metvaluesmax = log10(100);%fsc_options.conc_max_default);
  netgraph_concentrations(network,log10(my_c),[],1,gp); title('Optimised state: log10 Concentrations (mM)');
  
  figure(2); clf; 
  gp.arrowstyle='none'; gp.actstyle = 'fixed';
  gp.actvaluesmin = log10(0.01);
  gp.actvaluesmax = log10(100);
  gp.metvaluesmin = [];
  gp.metvaluesmax = [];
  netgraph_concentrations(fsc_options.network_CoHid,[],log10(my_A_forward),1,gp); title('Optimised state: log10 Reaction affinities (kJ/mol)');
  
  figure(3); clf; 
  gp.actnames  = gene_names;
  gp.actvaluesmin = log10(10^-6);
  gp.actvaluesmax = log10(10^-2);
  gp.metvaluesmin = [];
  gp.metvaluesmax = [];
  netgraph_concentrations(fsc_options.network_CoHid,[],log10(my_up),1,gp); title('Optimised state: log10 Enzyme levels (mM)');
  
  if fsc_options.print_graphics,
    if isfield(fsc_options,'psfile_dir'),
      cd(fsc_options.psfile_dir);
      print([ model_name '_' fsc_options.run_id '_' fsc_score '_opt_concentrations.eps'],'-f1','-depsc');
      print([ model_name '_' fsc_options.run_id '_' fsc_score '_opt_affinities.eps'],'-f2','-depsc');
      print([ model_name '_' fsc_options.run_id '_' fsc_score '_opt_enzymes.eps'],'-f3','-depsc');
    end
  end
  
end
