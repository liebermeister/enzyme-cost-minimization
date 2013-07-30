function cs_display_distribution_errorbar(fignum,mmatrix,names,ttitle,m_lower,m_upper,m_mean,options)

eval(default('options','struct'));
options_default = struct('flag_logarithmic',0,'shift',0);
options_default.colors = {'c','b','r'};
options = join_struct(options_default,options);

if isempty(m_lower), m_lower = nanmin(mmatrix(:)); end
if isempty(m_upper), m_upper = nanmax(mmatrix(:)); end

n_rel = size(mmatrix,1);
figure(fignum); subplot('position',[0.1 0.1 0.8 0.5]);
mmedian  = median(mmatrix')';
mquant5  = quantile(mmatrix',0.05)';
mquant95 = quantile(mmatrix',0.95)';
mmin     = min(mmatrix')';
mmax     = max(mmatrix')';
hold on
errorbar(options.shift+[1:n_rel],mmedian,mmedian-mmin,mmax-mmedian,'o','Color',options.colors{1});
errorbar(options.shift+[1:n_rel],mmedian,mmedian-mquant5,mquant95-mmedian,'o','Color',options.colors{2});
if ~iscell(m_mean),m_mean = {m_mean}; end
symbols= {'*','o'};
for it = 1:length(m_mean),
  m1  = m_mean{it};
  plot(find(isfinite(m1))+options.shift,m1(isfinite(m1)),'*','Color',options.colors{3},'Marker',symbols{it});
end
if options.flag_logarithmic, set(gca,'YScale','Log'); end 
axis([0, n_rel+1,m_lower,m_upper]);
a = axis;
my_xticklabel([1:n_rel]',a(4),names);
xlabel(ttitle);

