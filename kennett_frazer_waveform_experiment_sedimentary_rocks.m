close all
clear all

% Q estimation for waveform experiment from Kennett Frazer
vel_carbonate = 3200; % Adam et.al., 2009
vel_sandstone = 5400; % Subramaniyan + fluid substitution
den_carbonate = 2.368*1000; % Limestone + 20% porosity water
den_sandstone = 2.485 *1000; % Subramaniyan 
D = 52e-3; % Length of sample 52e-3
P1 = 0.5; % Proportion of plastic
P2 = 1-P1; % Proportion of steel
Q_carbonate =20; %1/Q = 0.05 carbonate Subramaniyan 2014 and Adam 2009
Q_sandstone = 50; %1/Q = 0.02 sandstone Subramaniyan 2015

M_period = [1 2 3 4 7 8 9 10 12 14 16 32 64 128 256]; % Periodicity
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
    
    vel = [vel_carbonate; vel_sandstone];
    
    den = [den_carbonate; den_sandstone];
    
    thick = [d1; d2];
    
    Q = [Q_carbonate; Q_sandstone];
    
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

%% Guide lines for the qplots

guide_invQ_visco = [0.256	0.041
0.359	0.041
0.754	0.040
0.863	0.041
0.903	0.041
0.964	0.043
1.116	0.047
1.306	0.051
1.606	0.054
1.865	0.056
2.235	0.057
2.699	0.056
2.740	0.055
2.778	0.052
3.013	0.049
3.307	0.047
4.061	0.046
5.063	0.045
6.540	0.044
12  0.0437
24 0.0436
48 0.0436
96 0.0436
];


guide_invQ_elas =[0.380	0.085
0.507	0.083
0.625	0.083
0.747	0.085
0.810	0.099
0.818	0.159
0.820	0.125
0.854	0.210
0.898	0.238
0.903	0.263
0.949	0.290
0.971	0.370
0.985	0.332
1.023	0.400
1.053	0.433
1.056	0.481
1.090	0.520
1.168	0.535
1.261	0.513
1.426	0.496
1.609	0.471
1.834	0.445
2.109	0.420
2.538	0.387
3.283	0.352
4.136	0.312
5.438	0.269
6.455	0.239
7.591	0.212
9.084	0.189
11.037	0.173
12.992	0.166
16.573	0.156
21.359	0.150
27.096	0.148
36.520	0.147
49.764	0.149
69.060	0.154
87.624	0.155
107.537	0.157
];

guide_invQ_elas(:,2) = guide_invQ_elas(:,2) ./100;

guide_invQ =[0.379	4.026
0.499	4.013
0.691	3.941
0.880	4.001
1.048	4.286
1.157	4.450
1.325	4.722
1.386	4.908
1.513	5.078
1.743	5.200
1.979	5.257
2.384	5.244
2.482	5.141
2.537	4.991
2.759	4.780
3.123	4.600
3.405	4.416
3.676	4.311
4.065	4.212
4.266	4.128
4.531	4.100
5.471	4.108
6.344	4.128
7.424	4.152
9.578	4.179
12.511	4.187
24 0.0421*100
48 0.0420*100
96 0.0420*100
];

guide_invQ(:,2) = guide_invQ(:,2) ./100;

%% Plotting Q curves for elastic, viscoelastic and intrinsic
fig1 = figure;
set(fig1, 'Units','inches', 'Position',[0 0 20 7],'PaperPositionMode','auto');
set(gca,...
'Units','normalized',... 
'FontUnits','points',... 
'FontWeight','normal',... 
'FontSize',16,... 
'FontName','Times')

rng(0,'twister');

invQ_disp_kf_elas = 1./Q_disp_kf_elas;
invQ_disp_kf_elas(invQ_disp_kf_elas<0) = (0.001.*rand(1));


subplot 131
semilogx(lambdaoverd_kf, 1./Q_disp_kf, 'ok');
hold on;
semilogx(guide_invQ_visco(:,1), guide_invQ_visco(:,2), '-k');
xlabel('\lambda/d');
ylabel('1/Q');
title('(i) Effective (Viscoelastic)');
box on; 
set(gca, 'Layer', 'top');
ylim([0.03 0.06]);
xlim([1e-1 1e3]);
set(gca, 'XTick', logspace(-1,3,5));
axis square;

subplot 132
% figure;
semilogx(lambdaoverd_kf_elas, invQ_disp_kf_elas, 'ok');
hold on;
semilogx(guide_invQ_elas(:,1), guide_invQ_elas(:,2), '-k');
semilogx(guide_invQ(:,1), guide_invQ(:,2), '-k');
xlabel('\lambda/d');
ylabel('1/Q');
title('(ii) Eff-scattering (elastic)');
box on; 
set(gca, 'Layer', 'top');
ylim([0 0.01]);
xlim([1e-1 1e3]);
set(gca, 'XTick', logspace(-1,3,5));
axis square;



subplot 133
% figure;
semilogx(lambdaoverd_kf, (1./Q_disp_kf)-invQ_disp_kf_elas, 'ok');
hold on;
semilogx(lambdaoverd_kf(end)+100, invq_eff_emt(end), '^k',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');
semilogx(guide_invQ(:,1), guide_invQ(:,2), '-k');
xlabel('\lambda/d');
ylabel('1/Q');
title('(iii) Effective - intrinsic');
box on; 
set(gca, 'Layer', 'top');
ylim([0.03 0.06]);
xlim([1e-1 1e3]);
set(gca, 'XTick', logspace(-1,3,5));
axis square;
legend('K-F Q^{-1}', 'Long wavelength limit Q^{-1}');

%% Plotting velocity dispersion curves from K-F
figure;
semilogx(lambdaoverd_kf, vel_disp_kf, 'o-k');
hold on;
semilogx(lambdaoverd_kf(end)+50, vel_eff_emt(end),'^k',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');
semilogx(lambdaoverd_kf(1)-0.2, real(vel_eff_rt(1)),'pk',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');
xlabel('\lambda/d');
ylabel('Velocity (m/sec)');
box on; 
set(gca, 'Layer', 'top');
ylim([3700 4200]);
xlim([1e-1 1e3]);
set(gca, 'XTick', logspace(-1,3,5));
legend('K-F velocities', 'Effective medium theory', 'Ray theory'); 
    