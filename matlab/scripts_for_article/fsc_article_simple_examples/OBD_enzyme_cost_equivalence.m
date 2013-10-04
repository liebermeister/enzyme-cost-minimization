% prove that obd and enzyme cost (for full saturation, all kcat the same)
% are not equivalent. plot both optimality criteria in the space of 
% thermo driving forces. constraint lines (edges of feasible polytope)
% may be parallel to the contour lines in different places for the different functions.

A_list = [0.01:0.1:5];

A1 = repmat(A_list,length(A_list),1);
A2 = A1';

min_A = [A1<=A2] .* A1 + [A1>A2] .* A2;

figure(1); contour(min_A); axis equal; axis square; xlabel('DF reaction 1'); ylabel('DF reaction 2');

enz_A = 1./[  1./[1-exp(-A1)] + 1./[1-exp(-A2)] ];

figure(2); contour(enz_A); axis equal; axis square; xlabel('DF reaction 1'); ylabel('DF reaction 2');

cd ~/projekte/obd/appendix_obd_mca/
print obd.eps -f1 -depsc   
print enz.eps -f2 -depsc
