% Velocity dispersion curve with different impedance contrasts

% close all; clear all;

vel_plastic = 2487; % 2487m/sec
vel_ratio = [1:0.1:3];
vel_steel = vel_ratio.*vel_plastic; % 5535m/sec
den_plastic = (1.741.*(vel_plastic./1000).^0.25).*1000; % Gardener's relation % 1.210 g/cc
den_steel = (1.741.*(vel_steel./1000).^0.25).*1000; % Gardener's relation % 7.900 g/cc
D = 52e-3; % Length of sample 52e-3
Q_plastic =10;
Q_steel = 20;

fdom = 200e3;

M = 1:256; % Periodicity

d = D./M; % Spatial period
n_layers = 2.*M; % Number of layers

vel_disp = zeros([length(vel_ratio),length(M)]);
Q_disp = zeros([length(vel_ratio),length(M)]);
lambdaoverd = zeros([length(vel_ratio),length(M)]);
vel_emt = zeros([length(vel_ratio),length(M)]);
vel_rt = zeros([length(vel_ratio),length(M)]);

for j = 1:length(vel_ratio)
    for i = 1:length(M)
        
        lyr = zeros(n_layers(i),3);
        
        vel = [vel_plastic; vel_steel(j)];
        
        den = [den_plastic; den_steel(j)];
        
        thick = [0.5.*d(i); 0.5.*d(i)];
        
        Q = [Q_plastic; Q_steel];
        
        lyr(:,1) = repmat(vel,[M(i),1]);
        lyr(:,2) = repmat(den,[M(i),1]);
        lyr(:,3) = repmat(thick,[M(i),1]);
        
        Q_layer = repmat(Q, [M(i),1]);
        
        lyr_kf = repmat(lyr,10,1);
        Q_layer_kf = repmat(Q_layer, 10,1);
        
        [freq_kf,vel_disp(j,i),~,~,Q_disp(j,i)] = kenfdispslowQ(lyr_kf,fdom,Q_layer_kf,2*pi*fdom);
        
        lambdaoverd(j,i) = (vel_disp(j,i)./freq_kf)./d(i);
        
        vel_emt(j,i) = veffemt(lyr,fdom,Q_layer,2*pi*fdom);
        vel_rt(j,i) = velrt_visco(lyr, fdom, Q_layer,2*pi*fdom);
        
    end
end


% vel_ratio_new = repmat(vel_ratio.',1,length(M));
% figure;
% contour(lambdaoverd(1:6,:), vel_ratio_new(1:6,:), vel_disp(1:6,:),30, 'LineWidth', 3);
% set(gca,'XScale', 'log');
% xlabel('\lambda/d');
% ylabel('V1/V2');
% c = colorbar('eastoutside');
% c.Label.String = 'Velocity (m/sec)';
% colormap('jet');
% xlim([1e-1 1e2]);


% figure;
% semilogx(lambdaoverd(1,:), (vel_disp(1,:)-real(vel_emt(1,:)))./(real(vel_rt(1,:))-real(vel_emt(1,:))), '-', 'LineWidth', 3);
% hold on;
% semilogx(lambdaoverd(6,:), (vel_disp(6,:)-real(vel_emt(6,:)))./(real(vel_rt(6,:))-real(vel_emt(6,:))), '-', 'LineWidth', 3);
% semilogx(lambdaoverd(11,:), (vel_disp(11,:)-real(vel_emt(11,:)))./(real(vel_rt(11,:))-real(vel_emt(11,:))), '-', 'LineWidth', 3);
% semilogx(lambdaoverd(21,:), (vel_disp(21,:)-real(vel_emt(21,:)))./(real(vel_rt(21,:))-real(vel_emt(21,:))), '-', 'LineWidth', 3);
% grid on; box on;
% xlabel('\lambda/d');
% ylabel('(V - V_{EMT})/(V_{RT} - V_{EMT})');
% legend(['Velocity ratio = ' num2str(vel_ratio(6))], ...
%     ['Velocity ratio = ' num2str(vel_ratio(11))],...
%     ['Velocity ratio = ' num2str(vel_ratio(21))]);
% ylim([-2 4]); xlim([1e-1 1e2]);
% 
% 
% figure;
% plot(vel_ratio, (real(vel_rt(:,1))-real(vel_emt(:,1)))./(real(vel_rt(:,1))), '-ok','MarkerSize',10, 'LineWidth',3);
% xlabel('Velocity ratio');
% ylabel('(V_{RT} - V_{EMT})/V_{RT}');
% title('(V_{RT} - V_{EMT})/V_{RT} vs Velocity ratio')
% grid on; box on;




%%
% Velocity dispersion curve with different Q contrasts

% close all; clear all;

vel_plastic = 2487; % 2487m/sec
% vel_ratio =2;
vel_steel = 2487; % 5535m/sec
den_plastic = 1.210.*1000; % Gardener's relation % 1.210 g/cc
den_steel =4.5.*1000; % Gardener's relation % 7.900 g/cc
D = 52e-3; % Length of sample 52e-3
Q_plastic =2.*pi;
Q_ratio = [1:1:100];
Q_steel = Q_plastic.*Q_ratio;

fdom = 200e3;

M = [1:1:256]; % Periodicity
M = 1;

d = D./M; % Spatial period
n_layers = 2.*M; % Number of layers

vel_disp = zeros([length(Q_ratio),length(M)]);
Q_disp = zeros([length(Q_ratio),length(M)]);
lambdaoverd = zeros([length(Q_ratio),length(M)]);

vel_disp_elas = zeros([length(Q_ratio),length(M)]);
Q_disp_elas = zeros([length(Q_ratio),length(M)]);

Qrt = zeros([length(Q_ratio),length(M)]);
invQemt = zeros([length(Q_ratio),length(M)]);

hmQ = zeros([length(Q_ratio),length(M)]); % Harmonic mean as effective Q


for j = 1:length(Q_ratio)
    for i = 1:length(M)
        
        lyr = zeros(n_layers(i),3);
        
        vel = [vel_plastic; vel_steel];
        
        den = [den_plastic; den_steel];
        
        thick = [0.5.*d(i); 0.5.*d(i)];
        
        Q = [Q_plastic; Q_steel(j)];
        
        lyr(:,1) = repmat(vel,[M(i),1]);
        lyr(:,2) = repmat(den,[M(i),1]);
        lyr(:,3) = repmat(thick,[M(i),1]);
        
        Q_layer = repmat(Q, [M(i),1]);
        
        lyr_kf = repmat(lyr,1,1);
        Q_layer_kf = repmat(Q_layer, 1,1);
        
        [freq_kf,vel_disp(j,i),~,~,Q_disp(j,i)] = kenfdispslowQ(lyr_kf,fdom,Q_layer_kf,2*pi*fdom);
        
        lambdaoverd(j,i) = (vel_disp(j,i)./freq_kf)./d(i);
        
        [~,Qrt(j,i)]=velrt_visco(lyr,fdom,Q_layer,2*pi*fdom);
        invQemt(j,i)=qeffemt(lyr,fdom,Q_layer,2*pi*fdom);
        
        frac = lyr(:,3)./sum(lyr(:,3));
        
        hmQ(j,i) = (sum(frac./Q_layer)).^-1;
        
        %         [freq_kf_elas,vel_disp_elas(j,i),~,~,Q_disp_elas(j,i)] = kenfdispslowQ(lyr_kf,fdom,Q_layer_kf.*1e10,2*pi*fdom);
        
    end
end


Q_ratio_new = repmat(Q_ratio.',1,length(M));
% figure;
% contour(lambdaoverd(1:6,:), 1./Q_ratio_new(1:6,:), 1./Q_disp(1:6,:),30, 'LineWidth', 3);
% set(gca,'XScale', 'log');
% xlabel('\lambda/d');
% ylabel('Q2/Q1');
% c = colorbar('eastoutside');
% c.Label.String = '1/Q';
% colormap('jet');
% xlim([1e-1 1e2]);


% figure;
% semilogx(lambdaoverd(1,:), (1./Q_disp(1,:)), '-', 'LineWidth', 3);
% hold on;
% semilogx(lambdaoverd(5,:), (1./Q_disp(5,:)), '-', 'LineWidth', 3);
% semilogx(lambdaoverd(10,:), (1./Q_disp(10,:)), '-', 'LineWidth', 3);
% grid on; box on;
% xlabel('\lambda/d');
% ylabel('1/Q');
% legend(['Q2/Q1 = ' num2str(Q_ratio(1))], ...
%     ['Q2/Q1 = ' num2str(Q_ratio(5))],...
%     ['Q2/Q1 = ' num2str(Q_ratio(10))]);
% xlim([1e-1 1e2]);
% ylim([0 0.4]);

% figure;
% semilogx(lambdaoverd(5,:), (1./Q_disp(1,:))-(1./Q_disp_elas(1,:)), '-', 'LineWidth', 3);
% hold on;
% semilogx(lambdaoverd(10,:), (1./Q_disp(5,:))-(1./Q_disp_elas(5,:)), '-', 'LineWidth', 3);
% semilogx(lambdaoverd(15,:), (1./Q_disp(15,:))-(1./Q_disp_elas(15,:)), '-', 'LineWidth', 3);
% grid on; box on;
% xlabel('\lambda/d');
% ylabel('1/Q');
% legend(['Q2/Q1 = ' num2str(Q_ratio(5))], ...
%     ['Q2/Q1 = ' num2str(Q_ratio(10))],...
%     ['Q2/Q1 = ' num2str(Q_ratio(15))]);
% xlim([1e-1 1e2]);
% % ylim([2500 4000]);


% figure;
% plot(Q_ratio, Qrt(:,1).*invQemt(:,1), '-ok', 'LineWidth', 3, 'MarkerSize',10);
% xlabel('Q2/Q1'); ylabel('Qrt/Qemt');
% grid on; box on;
% 
% figure;
% plot(Q_ratio, max(1./Q_disp,[],2)./(1./Q_plastic), '-ok','LineWidth', 3, 'MarkerSize',10);
% xlabel('Q2/Q1'); ylabel('(1/Q)_{max}/(1/Q1)');
% grid on; box on;
% 
% figure;
% plot(Q_ratio, invQemt./(1./Q_plastic), '-ok', 'LineWidth', 3, 'MarkerSize',10);
% hold on;
% plot(Q_ratio, (1./hmQ)./(1./Q_plastic), '-or', 'LineWidth', 3, 'MarkerSize',10);
% xlabel('Q2/Q1'); ylabel('(1/Q)_{eff}/(1/Q1)');
% grid on; box on;

%% Impedance contrast difference between Qeff and HM of Q

vel_plastic = 2487; % 2487m/sec
vel_ratio = [0.5:0.1:3];
vel_steel = vel_ratio.*vel_plastic; % 5535m/sec
den_plastic = (1.741.*(vel_plastic./1000).^0.25).*1000; % Gardener's relation % 1.210 g/cc
den_steel = (1.741.*(vel_steel./1000).^0.25).*1000; % Gardener's relation % 7.900 g/cc

impedance_ratio = (vel_plastic.*den_plastic)./(vel_steel.*den_steel);
D = 52; % Length of sample 52e-3
Q_plastic =100;
% Q_ratio = [0.1:0.1:2];
Q_ratio = logspace(-1,1,100);
Q_steel =Q_plastic.*Q_ratio;

fdom = 200e3;
M =1;
d = D./M; % Spatial period
n_layers = 2.*M; % Number of layers

invQemt = zeros([length(Q_ratio),length(vel_ratio)]);
% invQemt1 = zeros([length(Q_ratio),length(vel_ratio)]);

invQrt = zeros([length(Q_ratio),length(vel_ratio)]);
velemt_elastic =  zeros([length(Q_ratio),length(vel_ratio)]);
velemt_visco =zeros([length(Q_ratio),length(vel_ratio)]);

hmQ = zeros([length(Q_ratio),length(vel_ratio)]); % Harmonic mean as effective Q

for i = 1:length(Q_ratio)
for j = 1:length(vel_ratio)
    
    lyr = zeros(2,3);
    
    vel = [vel_plastic; vel_steel(j)];
    
    den = [den_plastic; den_steel(j)];
    
    thick = [0.5.*d; 0.5.*d];
    
    Q = [Q_plastic; Q_steel(i)];
    
    lyr(:,1) = vel;
    lyr(:,2) = den;
    lyr(:,3) = thick;
    
    Q_layer = Q;
    
%     invQemt(i,j)=qeffemt(lyr,fdom,Q_layer,2*pi*fdom);
    [invQemt(i,j),velemt_visco(i,j)] = viscoemt(lyr, Q_layer);
    [~,velemt_elastic(i,j)] = viscoemt(lyr, Q_layer.*1e100);
    
    frac = lyr(:,3)./sum(lyr(:,3));
    
    hmQ(i,j) = (sum(frac./Q_layer)).^-1;
    
    [~,Qrt] = velrt_visco(lyr,fdom,Q_layer,2*pi*fdom );
    invQrt(i,j) = 1./Qrt;
    
    
end
end

%%
contourlevels = [0.5:0.1:1.5 0.95 1.05];
figure;
[C,h] = contour(impedance_ratio, Q_ratio, (invQemt./(1./hmQ)),contourlevels, '-','Color', [0.83 0.82 0.78], 'LineWidth', 3, 'ShowText', 'on');
clabel(C,h, 'FontSize', 18);
xlabel('Impedance ratio');
ylabel('Q ratio');
title('Contours of Q^{-1}_{eff-intrinsic} / HM(Q^{-1})');
xlim([0.5 1.5]);
set(gca, 'YScale', 'log');
clabel(C,h, 'labelspacing', 200);

% contourlevels = 1:0.001:1.05;
% figure;
% [C,h] = contour(impedance_ratio, Q_ratio, (velemt_visco./velemt_elastic),contourlevels, '-','Color', [0.83 0.82 0.78], 'LineWidth', 3, 'ShowText', 'on');
% clabel(C,h, 'FontSize', 18);
% xlabel('Impedance ratio');
% ylabel('Q ratio');
% title('Contours of Q^{-1}_{eff-intrinsic} / HM(Q^{-1})');
% 
% contourlevels = [0.5:0.1:1, 2:1:9];
% figure;
% [C,h] = contour(impedance_ratio, Q_ratio, (invQemt)./(1./Q_plastic),contourlevels, '-','Color', [0.83 0.82 0.78], 'LineWidth', 3, 'ShowText', 'on');
% clabel(C,h, 'FontSize', 18);
% xlabel('Impedance ratio');
% ylabel('Q ratio');
% title('Contours of Q^{-1}_{eff-intrinsic} / HM(Q^{-1})');

% figure;
% plot(impedance_ratio, invQemt./(1./hmQ), '-k', 'LineWidth', 3);
% xlabel('Impedance ratio'); ylabel('(1/Q)_{eff}/HM(1/Q)');
% box on;
