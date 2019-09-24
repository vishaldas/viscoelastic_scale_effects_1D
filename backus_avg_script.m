% Running viscoelastic Backus Average toy problem
% Add srbtools to path before executing 

%% Welldata from GP262

clear all;
close all;
load('well2.mat');
Depth2 = well2.dept; %m
indexDepth = Depth2 >= 2100 & Depth2 <= 2296;
Depth2 = Depth2(indexDepth);
Rho2 = well2.rhob; %g/c3
Rho2 = Rho2*1000; %kg/m3 
Rho2 = Rho2(indexDepth);
Vp2 = well2.vp;  %km/s
Vp2 = Vp2*1000; %m/s
Vp2 = Vp2(indexDepth);
Vs2 = well2.vs;  %km/s
Vs2 = Vs2*1000; %m/s
Vs2 = Vs2(indexDepth);
GR2 = well2.gr; %GAPI
GR2 = GR2(indexDepth);
NPhi2 = well2.nphi;
NPhi2 = NPhi2(indexDepth); %Thermal Neutron Porosity (original Ratio Method) in Selected Lithology

%load corrected density and interpolate onto Depth2
load('Well2CorrectedDensity.mat');
load('Well2CorrectedDensityDepth.mat');
well2denscorr = well2denscorr*1000;
Density2 = interp1(well2denscorrdepth, well2denscorr, Depth2);

well2_sats = load('well_2_sats.txt');
well2_sats_depth = well2_sats(:,1); % in meters 

% Depth shift of 22 meters to match with the OWC and top of horizon
well2_sats_depth = well2_sats_depth +22;
well2_index3 = well2_sats_depth>=2100 & well2_sats_depth<=2296;
well2_sats_depth = well2_sats_depth(well2_index3);
well2_sats_deepsw = well2_sats(:,2); % in fractions
well2_sats_deepsw = well2_sats_deepsw(well2_index3);

load('denpor_brine.mat');
dencorr_phi = well2_dencorrpor_brine(well2_depth >= 2100 & well2_depth <= 2296);


%% Gassman's substitution

%Use Gassmann’s relations to do ?uid substitution in Well 2, ?rst from brine to
%oil, then from brine to gas in the upper 50 m of the reservoir sands. The ?uid
%properties are as follows: brine salinity 80,000 ppm, oil API 19, GOR 100 l/l,
%temperature 70 C, pore pressure 16 MPa, gas gravity 0.6.
BulkQuartz = 36.8e9;
ShearQuartz = 44e9;
BrineBulk = 2.8e9;
BrineDens = 1.090;
BrineBulk = 2.8; % GPa
OilDens = 0.780; 
OilGravity = 32; % API
GasGravity = 0.6;
GasIndexBrine = 0;
GasIndexOil = 0;
GOR = 64; 
method = 0;
BrineSalinty = 80000;
Temperature = 77.2;
SatOil = 1-well2_sats_deepsw;
SatGas = 0;

SatOil1 = 0;
SatGas1 = 1-well2_sats_deepsw;
% SatGas1(SatGas1>=0.6) = 0.6;
% SatGas1(SatGas1>=0.1) = 0.1; % Fizz gas


VpGassman = zeros(length(Depth2), 1);
VsGassman = zeros(length(Depth2), 1);
RhoGassman = zeros(length(Depth2), 1);


for i=1:length(Depth2)
    PorePressure = 20;
    % Initial saturation oil
    [Kreuss,rhoeff,~,~,Rhow,Kw,~,~,~,~,~,~,~]=flprop(method,...
    BrineSalinty,OilGravity,GasGravity,GOR,GasIndexBrine,GasIndexOil,...
    PorePressure,Temperature,SatOil(i),SatGas);
    
    % Final saturation gas
    [Kreuss1,rhoeff1,~,~,~,~,~,~,~,~,~,~,~]=flprop(method,...
    BrineSalinty,OilGravity,GasGravity,GOR,GasIndexBrine,GasIndexOil,...
    PorePressure,Temperature,SatOil1,SatGas1(i));


    [vpgas,vsgas,rogas,kgas]=gassmnv(Vp2(i),Vs2(i),Density2(i),...
        rhoeff*1000,Kreuss*1e9,rhoeff1*1000,Kreuss1*1e9,BulkQuartz,dencorr_phi(i));
    VpGassman(i) = vpgas;
    VsGassman(i) = vsgas;
    RhoGassman(i) = rogas;

end

% 
figure;

subplot 121
plot(Vp2,Depth2, '-k');
hold on;
plot(VpGassman, Depth2,'-b');
xlabel('Vp (m/s)');
ylabel('Depth (m)');
axis ij; 

subplot 122
plot(Density2,Depth2, '-k');
hold on;
plot(RhoGassman, Depth2,'-b');
xlabel('Density (kg/m3)');
ylabel('Depth (m)');
axis ij; 

%% Plotting the well data

well2_sats_deepsw = 1-SatGas1;
Vp2 = VpGassman;
Density2 = RhoGassman;
Vs2 = VsGassman;

figure;

subplot 151
plot(GR2, Depth2, '-k'); 
axis ij; 
xlabel('GR (gapi)'); 
ylabel('Depth (m)');
xlim([35 135]);

subplot 152
plot(Density2, Depth2, '-k'); 
axis ij; 
xlabel('Density (kg/m3)'); 

subplot 153
plot(Vp2, Depth2, '-k'); 
axis ij; 
xlabel('Vp (m/sec)'); 

subplot 154
plot(well2_sats_deepsw, Depth2, '-k'); 
axis ij; 
xlabel('Sw'); xlim([0 1]);
% 
% % Creating a synthetic invQ curve based on correlation with Sw log
% Qp = randwithcorr(well2_sats_deepsw,0.95,30,10);
% Qp(well2_sats_deepsw >=0.95) = 1e+4;
% Qp(Depth2>=2182) = 1e+4;
% Qp(Qp<=10) = 10;
% invQp = 1./Qp;

% Creating a invQ curve from Sw vs invQ curve
% Sw_lookup = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 0.965 0.98 1];
% Q_lookup = [300 200 100 50 40 30 25 20 10 15 20 50 100 300];
load('lookup.mat');
sw_interp = 0:0.01:1;
Q_lookup_interp = interp1(lookup(:,1),1./lookup(:,2),sw_interp);

Qp = spline(lookup(:,1),1./lookup(:,2),well2_sats_deepsw);
Qp(Qp<=9) = 9;
Qp(Depth2>=2182) = max(Qp); 

invQp = 1./Qp;

subplot 155
plot(invQp, well2_sats_depth, '-k'); 
axis ij; 
xlabel('Inv Qp'); 
% xlim([0 0.2]);
ylim([2100 2300]);

figure;
plot(sw_interp,1./Q_lookup_interp, '-k');
hold on;
% grid on; box on;
% plot(well2_sats_deepsw(Depth2<=2182),1./Qp(Depth2<=2182), 'ok');
xlabel('Sw'); ylabel('1/Q');
title('Sw vs 1/Q curve used for generating Q log');



%% Elastic Backus average 

% close all;

% window chosen considering seismic resolution of 4.5m
% [c_elast_bcks,rho_elast_bcks]=bkusrunavg(Depth2,...
%     Vp2,Vs2,Density2,30.*mean(diff(Depth2))); 
% 
% M_elast_bcks = c_elast_bcks(3,3,:);
% M_elast_bcks = M_elast_bcks(:);
% vp_elast_bcks = sqrt(M_elast_bcks./rho_elast_bcks);

[vp_elast_bcks2,vs_elast_bcks2,rho_elast_bcks2] = bkusrunavgnew(Depth2,Vp2,Vs2,Density2,30);

figure; 
subplot 121
plot(Density2, Depth2, '-k'); 
hold on;
plot(rho_elast_bcks2, Depth2, '--r'); 
% plot(rho_elast_bcks, Depth2, '--b'); 
axis ij; grid on; box on;
xlabel('Density (kg/m3)');
ylabel('Depth (m)');
legend('Original', 'Backus');
ylim([2130 2280]);

subplot 122
plot(Vp2, Depth2, '-k');
hold on;
plot(vp_elast_bcks2, Depth2, '--r'); 
% plot(vp_elast_bcks, Depth2, '--b'); 
axis ij; grid on; box on;
xlabel('Vp (m/sec)'); 
legend('Original', 'Backus');
ylim([2130 2280]);

%% Viscoelastic Backus average 


% window chosen considering seismic resolution of (10-15m)
% Sonic logging frequency at which the logs are specified is taken as 20KHz
% Upscaled to seismic frequency of 30 Hz 
% Viscoelastic model taken is Nearly constant Q model
[invqeff,vp_visco_bcks,rho_visco_bcks]=bkusrunavgvisco(Depth2,...
    Vp2,invQp,Density2,30,2.*pi.*30,2.*pi.*20000);

figure; 
subplot 131
plot(Density2, Depth2, '-k'); 
hold on;
plot(rho_visco_bcks, Depth2, '--r'); 
axis ij; grid on; box on;
xlabel('Density (kg/m3)');
ylabel('Depth (m)');
legend('Original', 'Backus');
ylim([2130 2280]);

subplot 132
plot(Vp2, Depth2, '-k');
hold on;
plot(vp_visco_bcks, Depth2, '--r'); 
axis ij; grid on; box on;
xlabel('Vp (m/sec)'); 
legend('Original', 'Backus');
ylim([2130 2280]);

subplot 133
plot(invQp, Depth2, '-k');
hold on;
plot(invqeff, Depth2, '--r'); 
axis ij; grid on; box on;
xlabel('invQp'); 
legend('Original', 'Backus');
ylim([2130 2280]);

%% Including only depths 2130 to 2280

Depth2_upscaled = Depth2(Depth2>=2130 & Depth2<=2280);
vp_elast_bcks2 = vp_elast_bcks2(Depth2>=2130 & Depth2<=2280);
vs_elast_bcks2 = vs_elast_bcks2(Depth2>=2130 & Depth2<=2280);
vp_visco_bcks = vp_visco_bcks(Depth2>=2130 & Depth2<=2280);
rho_elast_bcks2 = rho_elast_bcks2(Depth2>=2130 & Depth2<=2280);
rho_visco_bcks = rho_visco_bcks(Depth2>=2130 & Depth2<=2280);
invqeff = invqeff(Depth2>=2130 & Depth2<=2280);

sw_depth2 = well2_sats_deepsw;
well2_sats_deepsw = well2_sats_deepsw(Depth2>=2130 & Depth2<=2280);

%% Comparison between elastic and viscoelastic Backus average

figure;
subplot 121
plot(vp_elast_bcks2, Depth2_upscaled, '-k');
hold on
plot(vp_visco_bcks, Depth2_upscaled, '--r');
plot(linspace(1000,5000,100), linspace(2142,2142,100), '-b');
plot(linspace(1000,5000,100), linspace(2182,2182,100), '-b');
axis ij; grid on; box on;
xlabel('Vp (m/sec)'); 
ylabel('Depth');
legend('Elastic Backus', 'Viscoelastic Backus');
ylim([2130 2280]);
axis tight;

subplot 122
plot(vp_elast_bcks2.*rho_elast_bcks2, Depth2_upscaled, '-k');
hold on
plot(vp_visco_bcks.*rho_visco_bcks, Depth2_upscaled, '--r');
plot(linspace(1e6,1e7,100), linspace(2142,2142,100), '-b');
plot(linspace(1e6,1e7,100), linspace(2182,2182,100), '-b');
axis ij; 
xlabel('Ip (m/sec kg/m3)'); 
ylabel('Depth');
legend('Elastic Backus', 'Viscoelastic Backus');
axis tight;


%% Generating synthetic seismogram
% Source wavelet
fdom = 30; % Dominant frequency taken same as the upscaled log
dt_nyquist = 1/(2*fdom);
dt = 0.0351*dt_nyquist;
N_sample = round(dt_nyquist./dt);

% Depth to time
tt = 2*cumsum([diff(Depth2_upscaled); 0.1523]./vp_elast_bcks2');
tt = tt+2*(Depth2_upscaled(1)./vp_elast_bcks2(1));
ttinterp = min(tt):dt:max(tt);
vp_elast_bcks_time = interp1(tt,vp_elast_bcks2,ttinterp);
vs_elast_bcks_time = interp1(tt,vs_elast_bcks2,ttinterp);
vp_visco_bcks_time = interp1(tt,vp_visco_bcks,ttinterp);
rho_elast_bcks_time = interp1(tt,rho_elast_bcks2,ttinterp);
rho_visco_bcks_time = interp1(tt,rho_visco_bcks,ttinterp);
invqeff_time = interp1(tt,invqeff,ttinterp);
well2_sats_deepsw_time = interp1(tt,well2_sats_deepsw,ttinterp);

time = 0:dt:max(tt)*10;
wvlt = sourcewvlt;
wvlt(end:length(time)) = 0;

layer = zeros(length(vp_elast_bcks2),3);
layer(:,1) = vp_elast_bcks2';
layer(:,2) = rho_elast_bcks2';
layer(:,3) = [Depth2_upscaled(1) (Depth2_upscaled(2) ...
    - Depth2_upscaled(1)).*ones(1,length(vp_elast_bcks2)-1)];

[wz_elastic,~,~] = kennet(layer,wvlt,dt,2,1,-1); 

layer_visco = zeros(length(vp_visco_bcks),3);
layer_visco(:,1) = vp_visco_bcks';
layer_visco(:,2) = rho_visco_bcks';
layer_visco(:,3) = [Depth2_upscaled(1) (Depth2_upscaled(2) ...
    - Depth2_upscaled(1)).*ones(1,length(vp_visco_bcks)-1)];

[wz_visco,~,~] = kennettQ2(layer_visco,wvlt,dt,2,1,-1,(1./invqeff)',2.*pi.*fdom); 

%%

figure;

subplot 131
plot(well2_sats_deepsw_time, ttinterp, '-k');
hold on
plot(linspace(0,1,100), linspace(1.867,1.867,100), '-b');
plot(linspace(0,1,100),  linspace(1.893,1.893,100), '-b');
axis ij; 

subplot 132
plot(wz_elastic(time>=min(ttinterp) & time<=max(ttinterp)),...
    [time(time>=min(ttinterp) & time<=max(ttinterp))], '-k');
axis ij;
hold on;
plot([wz_visco(time>=min(ttinterp) & time<=max(ttinterp))],...
    [time(time>=min(ttinterp) & time<=max(ttinterp))], '-r');
legend('Elastic', 'Viscoelastic');
plot(linspace(-0.2,0.2,100), linspace(1.867,1.867,100), '-b');
plot(linspace(-0.2,0.2,100), linspace(1.893,1.893,100), '-b');


subplot 133
plot(vp_elast_bcks_time.*rho_elast_bcks_time, ttinterp, '-k');
hold on
plot(vp_visco_bcks_time.*rho_visco_bcks_time, ttinterp, '-r');
plot(linspace(1000,5000,100), linspace(1.867,1.867,100), '-b');
plot(linspace(1000,5000,100), linspace(1.893,1.893,100), '-b');
% [w,dt_ricker] = ricker(50);
% wt = 0:dt_ricker:30*dt_ricker;
% wti = 0:dt:30*dt_ricker;
% wi = interp1(wt,w,wti);
% r0 = [0 diff(log(vp_elast_bcks_time.*rho_elast_bcks_time))];
% ss = conv(r0,wi,'same');
% plot(ss, ttinterp, '--r');
plot(linspace(0,1,100), linspace(1.846,1.846,100), '-b');
plot(linspace(0,1,100),  linspace(1.876,1.876,100), '-b');
axis ij;

%% Plotting seismic
figure;
subplot 131
seisplot(wz_elastic(time>=min(ttinterp) & time<=max(ttinterp)),...
    min(ttinterp),dt);
ylim([1.85 1.95]);
ylabel('Time(sec)');
title('Elastic');
hold on;
plot(linspace(-25,25,100), linspace(1.867,1.867,100), '-r');
% plot(linspace(0,25,100),  linspace(1.893,1.893,100), '-r');

subplot 132
seisplot(wz_visco(time>=min(ttinterp) & time<=max(ttinterp)),...
    min(ttinterp)-0.012,dt);
ylim([1.85 1.95]);
set(gca, 'YTick', []);
ylabel('');
title('Visco-elastic');
hold on;
plot(linspace(-25,25,100), linspace(1.867,1.867,100), '-r');
% plot(linspace(0,25,100), linspace(1.893,1.893,100), '-r');


subplot 133
time_visco = time-0.012;
seisplot(wz_elastic(time>=min(ttinterp) & time<=max(ttinterp))...
    -wz_visco(time_visco>=min(ttinterp) & time_visco<=max(ttinterp)),...
    min(ttinterp),dt);
set(gca, 'YTick', []);
ylim([1.85 1.95]);
ylabel('');
hold on;
plot(linspace(-25,25,100), linspace(1.867,1.867,100), '-r');
title('Difference');


figure;
seisplot(wz_visco(time>=min(ttinterp) & time<=max(ttinterp)),...
    min(ttinterp)-0.012,dt);
ylim([1.85 1.95]);
ylabel('Time(sec)');
title('Normal-incidence seismogram');
hold on;
plot(linspace(-25,25,100), linspace(1.867,1.867,100), '-r');
% plot(linspace(0,25,100), linspace(1.893,1.893,100), '-r');

