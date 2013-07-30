% ------------------------------------------------------------------
% FSC figures for second and third metabolite in a chain c0 -> a -> b -> c3
% everything is computed explicitly
% ------------------------------------------------------------------

clear

conc_list = 10.^[-0.2:0.01:1];
nc     = length(conc_list);
kplus  = 1; 
kminus = 1;
c0     = 8;
c3     = 0.8;
v      = 1;

for it1 = 1:nc,
  for it2 = 1:nc
    a(it1,it2) = conc_list(it1);
    b(it1,it2) = conc_list(it2);
  end
end

w1_ma  = kplus * c0 - kminus * a;
w2_ma  = kplus * a  - kminus * b;
w3_ma  = kplus * b  - kminus * c3;

w1_cmr = w1_ma ./ [1 + c0 + a];
w2_cmr = w2_ma ./ [1 + a + b];
w3_cmr = w3_ma ./ [1 + b + c3];

u1 = 1./w1_cmr;
u2 = 1./w2_cmr;
u3 = 1./w3_cmr;
u1(find([w1_ma<0])) = nan;
u2(find([w2_ma<0])) = nan;
u3(find([w3_ma<0])) = nan;
u = u1 + u2 + u3;
u(find([w1_ma<0]+[w2_ma<0]+[w3_ma<0])) = nan;

figure(1); im(-1+log10(flipud(u1'))); axis square; xlabel('log a'); ylabel('log b'); title('Cost reaction 1');
figure(2); im(-1+log10(flipud(u2'))); axis square; xlabel('log a'); ylabel('log b'); title('Cost reaction 2');
figure(3); im(-1+log10(flipud(u3'))); axis square; xlabel('log a'); ylabel('log b'); title('Cost reaction 3');
figure(4); im(-1.5+log10(flipud(u' )),[-0.5 -0.2]); axis square; xlabel('log a'); ylabel('log b'); title('Total cost');

A1 = log([kplus/kminus]) - log(a./c0);
A2 = log([kplus/kminus]) - log(b./a);
A3 = log([kplus/kminus]) - log(c3./b);

% bei MTDF: loesung ist die gleiche egal welche monotone funktion angewendet wird!

mtdf = min(min(A1,A2),A3);
mtdf(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

sitdf = 1./A1 + 1./A2 + 1./A3;
sitdf(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

mitdf = max(max(1./A1, 1./A2), 1./A3);
mitdf(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

metdf = 1./(1-exp(-2*A1)) + 1./(1-exp(-2*A2)) + 1./(1-exp(-2*A3));
metdf(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

figure(5); im(1+log10(flipud(mtdf'))); axis square; xlabel('log a'); ylabel('log b'); title('MTDF');
figure(6); im(-1+log10(flipud(sitdf'))); axis square; xlabel('log a'); ylabel('log b'); title('Sum of inverse TDF');
figure(7); im(-1+log10(flipud(mitdf'))); axis square; xlabel('log a'); ylabel('log b'); title('Maximal inverse TDF');
figure(8); im(-1+log10(flipud(metdf')),[-0.5 ,-0.3]); axis square; xlabel('log a'); ylabel('log b'); title('Maximal predicted enzyme demand for full saturation');

% -----------------------------------------------------------------------------------

% cd /home/wolfram/projekte/flux_specific_cost/ps-files
% 
% print example_four_chain_u1.eps -f1 -depsc
% print example_four_chain_u2.eps -f2 -depsc
% print example_four_chain_u3.eps -f3 -depsc
% print example_four_chain_utot.eps -f4 -depsc
% print example_four_chain_mtdf.eps -f5 -depsc
% print example_four_chain_sidf.eps -f6 -depsc
% print example_four_chain_mitdf.eps -f7 -depsc
