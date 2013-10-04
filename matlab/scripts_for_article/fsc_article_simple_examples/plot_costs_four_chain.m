% ------------------------------------------------------------------
% FSC figures for second and third metabolite in a chain c0 -> a -> b -> c3
% everything is computed explicitly
% ------------------------------------------------------------------

clear

dd = 0.01;

enzyme_cost_colors = flipud([[1:-dd:0]' * [1 0.1 0] + [0:dd:1]' * [0.5 .5 1]; ...
                      [1:-dd:0]' * [0.5 0.5 1] + [0:dd:1]' * [0.95 0.95 1]]);

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

figure(1); set(gca,'Fontsize',24); 
im(-1+log10(flipud(u1')),[-1,2]); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Enzyme demand in reaction 1');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

figure(2); set(gca,'Fontsize',24); im(-1+log10(flipud(u2')),[-1,2]); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Enzyme demand in reaction 2');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

figure(3); set(gca,'Fontsize',24); im(-1+log10(flipud(u3')),[-1,2]); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Enzyme demand in reaction 3');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

figure(4); set(gca,'Fontsize',24); im(-1+log10(flipud(u' )),[-1,2]); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Total enzyme demand');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

% ----------------------------------------------------------------
% 3d graphics

figure(11);set(gca,'Fontsize',24);  axis equal
surf(a,b,log10(u1),'EdgeColor','none'); colormap(enzyme_cost_colors);
xlabel('log [A]'); ylabel('log [B]'); zlabel('log10 Enzyme demand');

figure(12);set(gca,'Fontsize',24);  axis equal 
surf(a,b,log10(u2),'EdgeColor','none');colormap(enzyme_cost_colors);
xlabel('log [A]'); ylabel('log [B]'); zlabel('log10 Enzyme demand');

figure(13);set(gca,'Fontsize',24);   axis equal
surf(a,b,log10(u3),'EdgeColor','none');colormap(enzyme_cost_colors);
xlabel('log [A]'); ylabel('log [B]'); zlabel('log10 Enzyme demand');

figure(14);set(gca,'Fontsize',24);   axis equal
surf(a,b,log10(u),'EdgeColor','none');colormap(enzyme_cost_colors);
xlabel('log [A]'); ylabel('log [B]'); zlabel('log10 Enzyme demand');



% ----------------------------------------------------------------
% bei OBD: loesung ist die gleiche egal welche monotone funktion angewendet wird!

A1 = log([kplus/kminus]) - log(a./c0);
A2 = log([kplus/kminus]) - log(b./a);
A3 = log([kplus/kminus]) - log(c3./b);


obd = min(min(A1,A2),A3);
obd(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

sibd = 1./A1 + 1./A2 + 1./A3;
sibd(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

mibd = max(max(1./A1, 1./A2), 1./A3);
mibd(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

mebd = 1./(1-exp(-2*A1)) + 1./(1-exp(-2*A2)) + 1./(1-exp(-2*A3));
mebd(find([A1<=0]+[A2<=0]+[A3<=0]))= nan;

figure(5); set(gca,'Fontsize',24); im(1+log10(flipud(obd'))); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title(sprintf('Minimal driving force'));
colormap([0.7 0.7 0.7; flipud(enzyme_cost_colors)]);

figure(6); set(gca,'Fontsize',24); im(-1+log10(flipud(sibd'))); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Total inverse driving force');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

figure(7); set(gca,'Fontsize',24); im(-1+log10(flipud(mibd'))); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Maximal inverse driving force');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

figure(8); set(gca,'Fontsize',24); im(-1+log10(flipud(mebd')),[-1,2]); axis square; 
xlabel('log [A]'); ylabel('log [B]'); title('Total enzyme demand (full saturation)');
colormap([0.7 0.7 0.7; enzyme_cost_colors]);

% -----------------------------------------------------------------------------------

cd /home/wolfram/projekte/flux_specific_enzyme_cost/ps-files/example_four_chain/

print example_four_chain_u1.eps -f1 -depsc
print example_four_chain_u2.eps -f2 -depsc
print example_four_chain_u3.eps -f3 -depsc
print example_four_chain_utot.eps -f4 -depsc

print example_four_chain_u1_3d.eps -f11 -depsc
print example_four_chain_u2_3d.eps -f12 -depsc
print example_four_chain_u3_3d.eps -f13 -depsc
print example_four_chain_utot_3d.eps -f14 -depsc

print example_four_chain_obd.eps -f5 -depsc
print example_four_chain_sidf.eps -f6 -depsc
print example_four_chain_mibd.eps -f7 -depsc
print example_four_chain_med.eps -f8 -depsc
