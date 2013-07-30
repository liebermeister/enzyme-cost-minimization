% ----------------------------------------------------------------------------------
% Compute and display enzyme costs of a simple mass action kinetics
% ----------------------------------------------------------------------------------

clear

conc_list = 10.^[-1:0.01:1];
nc = length(conc_list);

kplus = 1; 
kminus = 1;
keq = kplus/kminus;
v = 1;

for it1 = 1:nc,
  for it2 = 1:nc
    a(it1,it2) = conc_list(it1);
    b(it1,it2) = conc_list(it2);
  end
end

w_ma  = kplus * a - kminus * b;
w_cmr = w_ma ./ [1 + a + b];

w_ma_sinh  = - 2 * sqrt(kplus*kminus* [a .* b] ) .* sinh([log(keq) + log(b) - log(a)]/2);
w_cmr_sinh = w_ma_sinh ./ [1+a+b];

u_ma    = v./w_ma;
u_ma_2  = v.* [1 + 1./a] ./[1-exp([log(keq) + log(b) - log(a)])];
u_cmr   = v./w_cmr;
u_cmr_2 = v .* [1 + a + b]./a ./[1-exp([log(keq) + log(b) - log(a)])]; % exact
u_cmr_3 = v .* [1 + exp( - log(a)) ] ./[1 - exp([log(keq) + log(b) - log(a)])];     % approximated

t1 = [1 + exp( - log(a)) ];
t2 = 1./[1 - exp([log(keq) + log(b) - log(a)])];

ind_infeasible = find(b./a> keq);
w_cmr(ind_infeasible) = nan;
w_ma(ind_infeasible) = nan;
u_ma(ind_infeasible) = nan;
u_cmr(ind_infeasible) = nan;
u_ma_2(ind_infeasible) = nan;
u_cmr_2(ind_infeasible) = nan;
u_cmr_3(ind_infeasible) = nan;
t2(ind_infeasible) = nan;

% term_thermo = -sinh([log(keq) + log(b) - log(a)]/2); term_thermo(ind_infeasible) = nan;
% term_1_ma  = sqrt(a.*b);
% term_1_cmr = sqrt(a.*b)./(1+a+b);
%figure(1); im(log10(flipud(term_thermo')),[]); xlabel('a'); ylabel('b'); title('Thermoterm');
%figure(2); im(log10(flipud(term_1_ma')),[]); xlabel('a'); ylabel('b'); title('Prefactor term (ma)');
%figure(3); im(log10(flipud(term_1_cmr')),[]); xlabel('a'); ylabel('b'); title('Prefactor term (cmr)');
%figure(4); im(log10(flipud(u')),[]); xlabel('a'); ylabel('b'); title('Log10 enzyme cost (ma)');
%figure(6); im(log10(flipud(1./u')),[]); xlabel('a'); ylabel('log b'); title('Log10 rate (ma)');
%figure(7); im(log10(flipud((term_1_cmr.*term_thermo)')),[]); xlabel('a'); ylabel('log b'); title('Log10 rate (cmr)');

ca
figure(1); im(log10(flipud(u_cmr')),[]); xlabel('log a'); ylabel('log b');   title('Log10 enzyme cost  (cmr)'); colorbar
figure(2); im(log10(flipud(u_cmr_2')),[]); xlabel('log a'); ylabel('log b'); title('Log10 enzyme cost (cmr)');colorbar
figure(3); im(log10(flipud(u_cmr_3')),[]); xlabel('log a'); ylabel('log b'); title('Log10 enzyme cost, approximated (cmr)');colorbar
figure(12); im(log10(flipud(t1')),[]); xlabel('log a'); ylabel('log b'); title('Term 1');colorbar
figure(13); im(log10(flipud(t2')),[]); xlabel('log a'); ylabel('log b'); title('Term 2'); colorbar
figure(14); im(log10(flipud([t1.*t2]')),[]); xlabel('log a'); ylabel('log b'); title('Term 1*2'); colorbar

cd /home/wolfram/projekte/cba/optimal_concentration_profiles/ps-files
print example_uni_uni_reaction.eps -f2 -depsc

% ------------------------------------------------------------------

model_version = 'version_2010-12-21';   % 2011-02-26 nur zu testzwecken
filenames     = ycm_filenames(model_version);
load(filenames.network_file); 

[Mplus, Mminus, Wplus, Wminus, nm, nr] = make_structure_matrices(network.N,network.regulation_matrix,find(network.external));

for it = 1:nr,
  my_ind    = find(abs(Mplus(it,:))+abs(Mminus(it,:)));
  my_N      = network. N(my_ind,it)';
  my_Mplus  = Mplus(it,my_ind);
  B = 0.5;
  %B = 0.00001;
  %B = 0.9999;
  MM        = [my_Mplus' * my_Mplus] + B/[1-B]^2 * [my_N'*my_N] - B/[1-B] * [my_Mplus'*my_N+my_N'*my_Mplus];
  mmm = min(eig(MM));
  mmm(mmm>-10^-5) = 0;
 log10(svd(MM))   
 if mmm<0, display('buh!'); end
 pause
end
