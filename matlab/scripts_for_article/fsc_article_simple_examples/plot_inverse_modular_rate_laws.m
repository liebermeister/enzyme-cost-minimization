% -------------------------------------------------------------
% plot the inverse modular rate laws 
% -------------------------------------------------------------

clear

N = [-1 1]';
network = network_construct(N);

a_list = 10.^[-2.5:0.1:2.5];
b_list = 10.^[-2.5:0.1:2.5];

% -------------------------------------------------------------

network.kinetics = set_ms_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b]',network);
  end
end
v(v<=1/10) = nan;

figure(1);  colormap([ 0.8 0.9 1]);  set(gca,'FontSize',24);
h = surf(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v,'EdgeColor','blue');
hold on; plot3(a_list,a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.01 1 100],'YTick',[0.01 1 100],'XMinorGrid','off','YMinorGrid','off');
axis tight; %axis([0.01 100 0.01 100 -1 1]);
xlabel('a');ylabel('b');zlabel('v');
set(gca,'CameraPosition', [-6596.71 -3680.28 25.34164]);

% -------------------------------------------------------------

network.kinetics = set_cs_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b]',network);
  end
end
v(v<=1/10) = nan;

figure(2); colormap([ 0.8 0.9 1]); set(gca,'FontSize',24);
h = surf(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v,'EdgeColor','blue');
hold on; plot3(a_list,a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.01 1 100],'YTick',[0.01 1 100],'XMinorGrid','off','YMinorGrid','off'); 
axis tight; % axis([0.01 100 0.01 100 -1 1]);
xlabel('a');ylabel('b');zlabel('v');
set(gca,'CameraPosition', [-6596.71 -3680.28 25.34164]);

% -------------------------------------------------------------

N = [-1 1]';
network = network_construct(N);

a_list = 10.^[-1.5:0.1:1.5];
b_list = 10.^[-1.5:0.1:1.5];

network.kinetics = set_fd_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b]',network);
  end
end
v(v<=1/10) = nan;

figure(3); clf; colormap([ 0.8 0.9 1]); set(gca,'FontSize',24);
h = surf(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v,'EdgeColor','blue');
hold on; plot3(a_list,a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.1 1 10],'YTick',[0.1 1 10],'XMinorGrid','off','YMinorGrid','off'); 
axis tight; % axis([0.01 100 0.01 100 -1 1]);
xlabel('a');ylabel('b');zlabel('v');
set(gca,'CameraPosition', [-659 -368 302]);

% -------------------------------------------------------------

network.kinetics = set_rp_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b]',network);
  end
end
v(v<=1/10) = nan;

figure(4); colormap([ 0.8 0.9 1]); set(gca,'FontSize',24);
h = surf(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v,'EdgeColor','blue');
hold on; plot3(a_list,a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.1 1 10],'YTick',[0.1 1 10],'XMinorGrid','off','YMinorGrid','off'); 
axis tight; % axis([0.01 100 0.01 100 -1 1]);
xlabel('a');ylabel('b'); zlabel('v');
set(gca,'CameraPosition', [-659 -368 302]);

% -------------------------------------------------------------

N = [-1 -1 1]';
network = network_construct(N);

a_list = 10.^[-3:0.1:3];
b_list = 10.^[-3:0.1:3];

network.kinetics = set_cs_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b 1]',network);
  end
end
v(v<=1/10) = nan;

figure(5);  colormap([0 0 1]); set(gca,'FontSize',24);
h = mesh(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v);
hold on; plot3(a_list,1./a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.01 1 100],'YTick',[0.01 1 100]); 
axis tight; % axis([0.01 100 0.01 100 -1 1]);
xlabel('a'); ylabel('b');zlabel('v');
set(h,'FaceColor',[.9 .9 1],'FaceLighting','Phong')


network.kinetics = set_ds_kinetics(network);

clear v
for ita= 1:length(a_list),
  a = a_list(ita);
  for itb= 1:length(b_list),
    b = b_list(itb);
    v(ita,itb) = network_velocities([a b 1]',network);
  end
end
v(v<=1/10) = nan;

figure(6);  colormap([0 0 1]); set(gca,'FontSize',24);
h = mesh(repmat(a_list,length(b_list),1)',repmat(b_list,length(a_list),1),1./v);
hold on; plot3(a_list,1./a_list,0*a_list,'b--','LineWidth',3); hold off
set(gca,'FontSize',24,'XScale','log','YScale','log','XTick',[0.01 1 100],'YTick',[0.01 1 100]); 
axis tight; % axis([0.01 100 0.01 100 -1 1]);
xlabel('a'); ylabel('b');zlabel('v');
set(h,'FaceColor',[.9 .9 1],'FaceLighting','Phong')


if 0,
  cd ~/projekte/pathway_modelling/flux_specific_cost/ps-files/inverse_rate_laws
  print inverse_SMkinetics.eps -f1 -depsc
  print inverse_CMkinetics.eps -f2 -depsc
  print inverse_PMkinetics.eps -f3 -depsc
  print inverse_FMkinetics.eps -f4 -depsc
end

