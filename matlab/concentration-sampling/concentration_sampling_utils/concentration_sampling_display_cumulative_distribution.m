function concentration_sampling_display_cumulative_distribution(fignum,mmatrix,names,ttitle,mmin,mmax,options)

eval(default('options','struct'));
options_default = struct('flag_logarithmic',0);
options = join_struct(options_default,options);

if isempty(mmin), mmin = min(min(mmatrix)); end 
if isempty(mmax), mmax = max(max(mmatrix)); end 

[n_rel,n_sample] = size(mmatrix);
figure(fignum); clf; 
subplot('Position',[0.1 0.1 0.55 0.8]); 
cm = ryb_colors(n_rel);
[mmedian,order] = sort(median(mmatrix'));
for it = 1:n_rel,
  plot(sort(mmatrix(order(it),:)),[1:n_sample]/n_sample,'-','Color',cm(it,:)); hold on;
end
hold off;  axis([mmin,mmax,0,1]);
if options.flag_logarithmic, 
  set(gca,'XScale','Log');
end

h = legend(names(order),'FontSize',10);
set(h,'Position', [0.7 0.1 0.25 0.8]);
title(ttitle);
