% ------------------------------------------------------------------------
% Pareto triangle in flux space and enyzme space: 
%   simple example chain with triple branch point
% ------------------------------------------------------------------------

cd ~/matlab/projects/pathway_modelling/flux_specific_cost/pareto
load pareto_branch_network

c_init = [10 1 0.1 0.1 0.1]';

network.kinetics.KM  = 2 * sign(network.kinetics.KM);% .* exp(randn(size(network.kinetics.KM)));
network.kinetics.Keq = [1 0.1 1 10]';


% -----------------------------------------
% create mesh

[nodes,triangles] = make_triangle_mesh(4);


% -----------------------------------------
% run flux-specific cost analysis 

n_nodes = size(nodes,2);

V = network_analyse(network);

clear c_fsc2 c_fsc3 u_fsc2 u_fsc3

options.show_graphics       = 0;
options.flag_given_kinetics = 1;
options.enzyme_cost_weights = ones(4,1);%[1 2 3 4]';
options.fsc_scores          = {'fsc2','fsc3'};
options.met_fix             = network.metabolites(find(network.external));
options.conc_fix            = c_init(find(network.external));

for it = 1:n_nodes,
  display(sprintf('Node %d/%d',it,n_nodes))
  v             = V * nodes(:,it);
  [c_opt,u_opt] = flux_specific_enzyme_cost(network,v,options);
  c_fsc2(:,it)  = c_opt.fsc2;
  c_fsc3(:,it)  = c_opt.fsc3;
  u_fsc2(:,it)  = u_opt.fsc2;
  u_fsc3(:,it)  = u_opt.fsc3;
end

u_tot_fsc2 = options.enzyme_cost_weights' * u_fsc2;
u_tot_fsc3 = options.enzyme_cost_weights' * u_fsc3;


% -----------------------------------------
% graphics

figure(11); plot([u_fsc2(1,:); u_fsc2(2:end,:)./nodes; c_fsc2(2,:)]');
legend({'u1','u2/v2','u3/v3','u4/v4','c'});

figure(10);
netgraph_concentrations(network,network.external,[],1);

figure(1); % elementary mode coefficients
trimesh(triangles', nodes(1,:), nodes(2,:), nodes(3,:),'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.2 0.2 0.2],'FaceLighting','flat');
title('Coefficients of elementary modes'); set(gca,'CameraPosition',[4.5, 1.9 1.75]); axis square

figure(2); % enzyme levels
% and internal concentration (according to FSC2 measure)
trimesh(triangles', u_fsc2(2,:), u_fsc2(3,:), u_fsc2(4,:),c_fsc2(2,:),'FaceVertexCData',c_fsc2(2,:)','EdgeColor',[1 1 1],'FaceColor','interp');
title('Enzyme levels and internal concentration (FSC2)'); set(gca,'CameraPosition',[4.5, 1.9 1.75]); axis square
colormap(ryb_colors); colorbar

figure(3); % enzyme levels
% and internal concentration (according to FSC3 measure)
trimesh(triangles', u_fsc3(2,:), u_fsc3(3,:), u_fsc3(4,:),c_fsc3(2,:),'FaceVertexCData',c_fsc3(2,:)','EdgeColor',[1 1 1],'FaceColor','interp');
title('Enzyme levels and internal concentration (FSC3)'); set(gca,'CameraPosition',[4.5, 1.9 1.75]); axis square
colormap(ryb_colors); colorbar

figure(12); % enzyme levels
% and internal concentration (according to FSC2 measure)
trimesh(triangles', u_fsc2(2,:), u_fsc2(3,:), u_fsc2(4,:),u_tot_fsc2,'FaceVertexCData',u_tot_fsc2','EdgeColor',[1 1 1],'FaceColor','interp');
title('Enzyme levels and total cost (FSC2)'); set(gca,'CameraPosition',[4.5, 1.9 1.75]); axis square
colormap(ryb_colors); colorbar

figure(13); % enzyme levels
% and internal concentration (according to FSC3 measure)
trimesh(triangles', u_fsc3(2,:), u_fsc3(3,:), u_fsc3(4,:),u_tot_fsc3,'FaceVertexCData',u_tot_fsc3','EdgeColor',[1 1 1],'FaceColor','interp');
title('Enzyme levels and total cost (FSC3)'); set(gca,'CameraPosition',[4.5, 1.9 1.75]); axis square
colormap(ryb_colors); colorbar

if 0,
  cd  /home/wolfram/projekte/pathway_modelling/flux_specific_cost/ps-files/pareto
  print pareto_branch_network.eps -f10 -depsc  
  print pareto_branch_efficiencies.eps -f11 -depsc
  print pareto_branch_elemenatary_modes.eps -f1 -depsc
  print pareto_branch_fsc2_concentration.eps -f2 -depsc
  print pareto_branch_fsc3_concentration.eps -f3 -depsc
  print pareto_branch_fsc2_cost.eps -f12 -depsc
  print pareto_branch_fsc3_cost.eps -f13 -depsc
end
