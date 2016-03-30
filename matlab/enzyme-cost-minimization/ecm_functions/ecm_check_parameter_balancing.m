function ecm_check_parameter_balancing(r, r_orig, network, show_graphics)

eval(default('show_graphics','1'));

threshold_mu = 5; % kJ/mol
threshold = 2;

mu0_change = r.mu0-r_orig.mu0; 
ind_change = find(isfinite(mu0_change) .* abs(mu0_change)>0);
if length(ind_change),
 display(sprintf('Changes of mu0 values (additive change > %f kK/mol shown',threshold_mu));
 print_matrix(mu0_change(ind_change), network.metabolites(ind_change))
end

log_Keq_change = log10(r.Keq./r_orig.Keq); 
ind_change = find(isfinite(log_Keq_change) .* abs(log_Keq_change) > log10(threshold));
if length(ind_change),
display(sprintf('Changes of Keq values (fold change > %f shown)',threshold));
print_matrix(10.^log_Keq_change(ind_change), network.actions(ind_change))
end

log_Kcatf_change = log10(r.Kcatf./r_orig.Kcatf); 
ind_change = find(isfinite(log_Kcatf_change) .* abs(log_Kcatf_change) > log10(threshold));
if length(ind_change),
display(sprintf('Changes of Kcatf values (fold change > %f shown)',threshold));
print_matrix(10.^log_Kcatf_change(ind_change), network.actions(ind_change))
end

log_Kcatr_change = log10(r.Kcatr./r_orig.Kcatr); 
ind_change = find(isfinite(log_Kcatr_change) .* abs(log_Kcatr_change) > log10(threshold));
if length(ind_change),
display(sprintf('Changes of Kcatr values (fold change > %f shown)',threshold));
print_matrix(10.^log_Kcatr_change(ind_change), network.actions(ind_change))
end

log_KM_change = log10(full(r.KM./r_orig.KM));
log_KM_change(~isfinite(log_KM_change)) = 0;
indices = find(abs(log_KM_change(:))>log10(threshold));
if length(indices), 
display(sprintf('Changes of KM values (fold change > %f shown)',threshold));
[~,order] = sort(abs(log_KM_change(indices)));
indices = indices(order(end:-1:1));
[ind_i,ind_j] = ind2sub(size(log_KM_change), indices);
for it = 1:length(indices),
 display(sprintf(' %s / %s: %f', network.actions{ind_i(it)}, network.metabolites{ind_j(it)}, 10.^log_KM_change(indices(it)))); 
end
end

log_KA_change = log10(full(r.KA./r_orig.KA));
log_KA_change(~isfinite(log_KA_change)) = 0;
indices = find(abs(log_KA_change(:))>log10(threshold));
[~,order] = sort(abs(log_KA_change(indices)));
indices = indices(order(end:-1:1));
if length(indices), 
display(sprintf('Changes of KA values (fold change > %f shown)',threshold));
[ind_i,ind_j] = ind2sub(size(log_KA_change), indices);
for it = 1:length(indices),
 display(sprintf(' %s / %s: %f', network.actions{ind_i(it)}, network.metabolites{ind_j(it)}, 10.^log_KA_change(indices(it)))); 
end
end

log_KI_change = log10(full(r.KI./r_orig.KI));
log_KI_change(~isfinite(log_KI_change)) = 0;
indices = find(abs(log_KI_change(:))>log10(threshold));
[~,order] = sort(abs(log_KI_change(indices)));
indices = indices(order(end:-1:1));
if length(indices), 
display(sprintf('Changes of KI values (fold change > %f shown)',threshold));
[ind_i,ind_j] = ind2sub(size(log_KI_change), indices);
for it = 1:length(indices),
 display(sprintf(' %s / %s: %f', network.actions{ind_i(it)}, network.metabolites{ind_j(it)}, 10.^log_KI_change(indices(it)))); 
end
end

if show_graphics,
  figure(1); 
  subplot(3,3,1); plot(r_orig.mu0,  r.mu0,'.'); title('mu0');axis tight
  subplot(3,3,2); plot(r_orig.Keq,  r.Keq,'.'); title('Keq'); set(gca, 'XScale','log','YScale','log');axis tight; 
  subplot(3,3,3); plot(r_orig.Kcatf,r.Kcatf,'.'); title('Kcatf'); set(gca, 'XScale','log','YScale','log');axis tight
  subplot(3,3,4); plot(r_orig.Kcatr,r.Kcatr,'.'); title('Kcatr'); set(gca, 'XScale','log','YScale','log');axis tight
  subplot(3,3,5); plot(r_orig.KM(:),r.KM(:),'.'); title('KM'); set(gca, 'XScale','log','YScale','log');axis tight
  subplot(3,3,6); plot(r_orig.KA(:),r.KA(:),'.'); title('KA'); set(gca, 'XScale','log','YScale','log');axis tight
  subplot(3,3,7); plot(r_orig.KI(:),r.KI(:),'.'); title('KI'); set(gca, 'XScale','log','YScale','log'); axis tight
  xlabel('Input data'); ylabel('Balanced');
end 