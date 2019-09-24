%% Q estimation for waveform experiment from Kennett Frazer

clear all;
% close all;

% cd('/Users/vdas2/GoogleDrive_stanford/Research/Codes/EffectiveQ');
addpath('/Users/vdas2/GoogleDrive_stanford/Research/Codes/srbtools/');


vel_plastic = 2487; % 2487m/sec
vel_steel = 5535; % 5535m/sec
den_plastic = 1.21*1000; % 1.210 g/cc
den_steel = 7.9 *1000; % 7.900 g/cc
D = 52e-3; % Length of sample 52e-3
P1 = 0.5; % Proportion of plastic
P2 = 1-P1; % Proportion of steel
Q_plastic =10; %10
Q_steel = 20; %20

M_period = [1 2 3 5 6 7 8 9 10 12 14 16 32 64 128 256]; % Periodicity

% M_period = [1 2 3 5 6 8 9 10 12 14 16 32 64 128 256]; % Periodicity

% M_period = 10;

fdom = 200e3;
fdom1 = 115e3; %115e3
% freq_kf = 10.^([0:0.1:10]);

vel_disp_kf = zeros(length(M_period),1);
Q_disp_kf = zeros(length(M_period),1);
lambdaoverd_kf = zeros(length(M_period),1);

vel_disp_kf_elas = zeros(length(M_period),1);
Q_disp_kf_elas = zeros(length(M_period),1);
lambdaoverd_kf_elas = zeros(length(M_period),1);

vel_eff_emt = zeros(length(M_period),1);
invq_eff_emt = zeros(length(M_period),1);

vel_eff_emt_elas = zeros(length(M_period),1);
invq_eff_emt_elas = zeros(length(M_period),1);

vel_eff_rt = zeros(length(M_period),1);
q_eff_rt = zeros(length(M_period),1);

vel_eff_rt_elas = zeros(length(M_period),1);
q_eff_rt_elas = zeros(length(M_period),1);

for i = 1:length(M_period)
    disp (['Iteration no ' num2str(i) '/' num2str(length(M_period))]);
    
    M = M_period(i);
    
    n_layers = 2*M; % Number of layers
    
    d = D/M; % Spatial period
    
    d1 = P1*d; %thickness of plastic
    d2= P2*d; %thickness of steel
    
    lyr = zeros(n_layers,3);
    
    vel = [vel_plastic; vel_steel];
    
    den = [den_plastic; den_steel];
    
    thick = [d1; d2];
    
    Q = [Q_plastic; Q_steel];
    
    lyr(:,1) = repmat(vel,[M,1]);
    lyr(:,2) = repmat(den,[M,1]);
    lyr(:,3) = repmat(thick,[M,1]);
    
    Q_layer = repmat(Q, [M,1]);
    
    lyr_kf = repmat(lyr,10,1);
    Q_layer_kf = repmat(Q_layer, 10,1);
    
    [fdom,vel_disp_kf(i),~,~,Q_disp_kf(i)] = kenfdispslowQ(lyr_kf,fdom,Q_layer_kf,2*pi*fdom1);
    lambdaoverd_kf(i) = (real(vel_disp_kf(i))./fdom)./d;
    
    vel_eff_emt(i) = veffemt(lyr_kf,fdom,Q_layer_kf,2*pi*fdom1);
    invq_eff_emt(i) = qeffemt(lyr_kf,fdom,Q_layer_kf,2*pi*fdom1);

    [vel_eff_rt(i), q_eff_rt(i)] = velrt_visco(lyr_kf,fdom,Q_layer_kf,2*pi*fdom1);
    
%     
%     [fdom,vel_disp_kf_elas(i),~,~,Q_disp_kf_elas(i)] = kenfdispslowQ(lyr_kf,fdom,Q_layer_kf.*1e10,2*pi*fdom);
%     lambdaoverd_kf_elas(i) = (real(vel_disp_kf_elas(i))./fdom)./d;
    
%     vel_eff_emt_elas(i) = veffemt(lyr_kf,fdom,Q_layer_kf.*1e10,2*pi*fdom);
%     invq_eff_emt_elas(i) = qeffemt(lyr_kf,fdom,Q_layer_kf.*1e10,2*pi*fdom);
    [vel_eff_rt_elas(i), q_eff_rt_elas(i)] = velrt_visco(lyr_kf,fdom,Q_layer_kf.*1e10,2*pi*fdom1);
    
    [fdom,vel_disp_kf_elas(i),~,~,Q_disp_kf_elas(i)] = kenfdispslow(lyr_kf,fdom);
    lambdaoverd_kf_elas(i) = (real(vel_disp_kf_elas(i))./fdom)./d;
    
    
    
    
%     [freq_kf,vel_disp_kf,~,~,Q_disp_kf] = kenfdispslowQ(lyr_kf,freq_kf,Q_layer_kf,2*pi*fdom);
%     vel_eff_emt = veffemt(lyr_kf,freq_kf,Q_layer_kf,2*pi*fdom);
%     invq_eff_emt = qeffemt(lyr_kf,freq_kf,Q_layer_kf,2*pi*fdom);
%     
%     vel_disp_kf(freq_kf <= 500) = vel_eff_emt(freq_kf <=500);
%     Q_disp_kf(freq_kf <= 500) = 1./(invq_eff_emt(freq_kf<=500));`
%     
%     lambdaoverd_kf = (real(vel_disp_kf)./freq_kf)./d;
%     
%     [freq_kf,vel_disp_kf_elas,~,~,Q_disp_kf_elas] = kenfdispslow(lyr_kf,freq_kf);
%     
%     lambdaoverd_kf_elas = (real(vel_disp_kf_elas)./freq_kf)./d;
    
end

%% Calculating full curve for K-F

freq_fkf = sort([1e-3 logspace(2,10,100) 200e3 115e3]);

[freq_fkf,veldisp1_fkf, ~, ~,Qdisp1_fkf]=kenfdispslowQ(lyr,freq_fkf, Q_layer,2*pi*fdom1);
q_emt1 = 1./(qeffemt(lyr,freq_fkf,Q_layer,2*pi*fdom));
lambdaoverd_fkf = real(veldisp1_fkf./freq_fkf)./d;

[freq_fkf,veldisp1_fkf_elas,~,~,Qdisp1_fkf_elas] = kenfdispslow(lyr,freq_fkf);
lambdaoverd_fkf_elas = real(veldisp1_fkf_elas./freq_fkf)./d;

%% Guide curve
guide_invQ = [
0.313	0.051
0.480	0.051
0.691	0.051
0.848	0.049
0.941	0.044
1.006	0.035
1.044	0.024
1.400	0.011
1.768	0.008
2.017	0.024
2.002	0.044
2.020	0.072
2.073	0.117
2.039	0.143
2.038	0.178
2.061	0.209
2.418	0.181
2.644	0.160
2.773	0.145
3.079	0.129
3.599	0.117
4.592	0.112
7.644	0.107
10.139	0.105
13.807	0.101
21.554	0.096
22.3    0.0966
45  0.1020
];

%% Plotting Q curves for elastic, viscoelastic and intrinsic
fig1 = figure;
set(fig1, 'Units','inches', 'Position',[0 0 20 7],'PaperPositionMode','auto');
set(gca,...
'Units','normalized',... 
'FontUnits','points',... 
'FontWeight','normal',... 
'FontSize',16,... 
'FontName','Times')



subplot 131
% semilogx(lambdaoverd_fkf, 1./Qdisp1_fkf, '--k');
semilogx(lambdaoverd_kf, 1./Q_disp_kf, 'ok');
hold on;
% Fitting a Lorentzian curve as a guide
[yprime2, P] = lorentzfit(lambdaoverd_kf,1./Q_disp_kf,[],[],'3c');
xx = logspace(log10(0.3),log10(50) , 100);
guide_Q_disp_kf = P(1)./((xx - P(2)).^2 + P(3)) + P(4);
semilogx(xx, guide_Q_disp_kf, '-k');
xlabel('\lambda/d');
ylabel('1/Q');
title('(i) Effective (Viscoelastic)');
box on; 
set(gca, 'Layer', 'top');
ylim([0 3.5]);
xlim([1e-1 1e2]);
set(gca, 'XTick', logspace(-1,2,4));
axis square;

subplot 132
% figure;
% semilogx(lambdaoverd_fkf_elas(1./Qdisp1_fkf_elas >= 0.004), 1./Qdisp1_fkf_elas(1./Qdisp1_fkf_elas >= 0.004), '--k');
semilogx(lambdaoverd_kf_elas, 1./Q_disp_kf_elas, 'ok');
hold on;
% Fitting a Lorentzian curve as a guide
invQ_disp_kf_elas = 1./Q_disp_kf_elas;
invQ_disp_kf_elas(invQ_disp_kf_elas<=0.1) = 0.1;
% lambdaoverd_kf_elas = lambdaoverd_kf_elas(invQ_disp_kf_elas>0);
[yprime2, P] = lorentzfit(lambdaoverd_kf_elas,invQ_disp_kf_elas,P,[],'3c');
xx = logspace(log10(0.3),log10(50) , 20);
guide_Q_disp_kf = P(1)./((xx - P(2)).^2 + P(3)) + P(4);
semilogx(xx(guide_Q_disp_kf<=3.5 & guide_Q_disp_kf>=0), ...
    guide_Q_disp_kf(guide_Q_disp_kf<=3.5 & guide_Q_disp_kf>=0), '-k');
xlabel('\lambda/d');
ylabel('1/Q');
title('(ii) Eff-scattering (elastic)');
box on; 
set(gca, 'Layer', 'top');
ylim([0 3.5]);
xlim([1e-1 1e2]);
set(gca, 'XTick', logspace(-1,2,4));
axis square;



subplot 133
% figure;

semilogx(lambdaoverd_kf, (1./Q_disp_kf)-(1./Q_disp_kf_elas), 'ok');
hold on;
semilogx(lambdaoverd_kf(end)+20, invq_eff_emt(end), '^k',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');

% Plotting a smooth curve through points
% % Selecting points manually to make a smooth curve 
% y = (1./Q_disp_kf)-(1./Q_disp_kf_elas);
% vec = [1:7 13:16];
% xx = lambdaoverd_kf(vec);
% yy = y(vec);
% xx = [xx; 0.95; 4]; yy = [yy; 0.03; 0.125];
% p = spline(xx,yy);
% xx1 = logspace(log10(0.3),log10(50) , 20);
% v = ppval(p, xx1);
% plot(xx1(v <= max(y)),v(v<=max(y)), '-k');

plot(guide_invQ(:,1), guide_invQ(:,2), '-k');

xlabel('\lambda/d');
ylabel('1/Q');
title('(iii) Effective - intrinsic');
box on; 
set(gca, 'Layer', 'top');
ylim([0 0.3]);
xlim([1e-1 1e2]);
set(gca, 'XTick', logspace(-1,2,4));
axis square;
legend('K-F Q^{-1}', 'Long wavelength limit Q^{-1}');
    