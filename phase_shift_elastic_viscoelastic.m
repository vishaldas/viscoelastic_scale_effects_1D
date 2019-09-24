clear all;
% close all;

addpath('/Users/vdas2/GoogleDrive_stanford/Research/Codes/srbtools/');
addpath(genpath('/Users/vdas2/GoogleDrive_stanford/Research/Codes/crewes/'));

% Three layered model with middle layer as viscoelastic
vel = [3000; 3000; 3000];
den = [2350; 2350; 2350];
thick = [500; 250; 500];
Q = [100; 10; 50];

ray_time = 2.*sum(thick./vel);

lyr =  [vel den thick];

fdom = 30;
dt_nyquist = 1/(2*fdom);
dt = 0.0351*dt_nyquist;
% [wvlt,time]= ricker_crewes(dt,fdom, 16);
% wvlt = wvlt(time >= 0).';
% time =  time(time >=0).';

wvlt = sourcewvlt();
time = 0:dt:7;
wvlt(end:length(time)) = 0;

[wz,pz,tf]=kennettQ2(lyr,wvlt,dt,2,1,-1,Q,2*pi*fdom);
[wz_elas,pz_elas,tf_elas]=kennet(lyr,wvlt,dt,2,1,-1);

% muting first few samples as they are the direct record
wz = wz(find(time >= 0.1,1)).*(time <= 0.1) + wz.*(time>0.1);
wz_elas = wz_elas(find(time >= 0.1,1)).*(time <= 0.1) + wz_elas.*(time>0.1);

figure;
subplot 131
stairs([0;vel], [0;cumsum(thick)], '-k', 'LineWidth', 3);
hold on;
stairs([0;den], [0;cumsum(thick)], '--k', 'LineWidth', 3);
% stairs([0;10000./Q], [0;cumsum(thick)], '-r', 'LineWidth', 3);
axis ij; box on;
set(gca, 'Layer', 'top'); % To bring axis on top
ylim([1,sum(thick)]);
xlim([min(den)-100 max(vel)+100]);
legend('Vp (m/sec)', 'Density (kg/m^3)');
title('Model', 'FontWeight', 'bold');
ylabel('Depth (m)');

subplot 132
stairs([0;1./Q], [0;cumsum(thick)], '-k', 'LineWidth', 3);
axis ij; box on;
set(gca, 'Layer', 'top'); % To bring axis on top
ylim([1,sum(thick)]);
% xlim([0 max(vel)+1000]);
title('Model', 'FontWeight', 'bold');
legend('Qp ^{-1}');
ylabel('Depth (m)');
xlabel('Qp ^{-1}');
set(gca,'Xtick',0:0.02:0.1)


subplot 133
plot(pz, time, '-k');
hold on;
plot(pz_elas, time, '--k');
ylim([0.35 0.5]);
xlim([-1 1.1]);
axis ij;
legend('Viscoelastic', 'Elastic');
title('Transmitted waveforms', 'FontWeight', 'bold');
ylabel('Time (sec)');
