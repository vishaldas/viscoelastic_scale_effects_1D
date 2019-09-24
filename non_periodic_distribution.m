%% Generating the random medium

close all;
clear all;

addpath('/Users/vdas2/GoogleDrive_stanford/Research/Codes/srbtools');

rng('default');

% ------------------------- Material properties ---------------------------
vel_plastic = 2487; % 2487m/sec
vel_steel = 5535; % 5535m/sec
velavg = (vel_plastic+vel_steel)/2;
den_plastic = 1.210*1000; % 1.210 g/cc
den_steel = 7.900 *1000; % 7.900 g/cc
thick_plastic = 540e-6; %540e-6
thick_steel = 490e-6;
Q_plastic = 10; % Quality factor for plastic Q=10
Q_steel = 20; % Quality factor for steel Q=20
Q_avg = (Q_plastic+Q_steel)/2;

den_avg = (den_plastic+den_steel)/2;

% % ---------------------- Exponential sequence ---------------------------- 
% n_r = 100; % Number of realizations
% n = 200; % number of layers
% corr = 7; % Correlation length of exponential function
% vel = zeros(n,n_r);
% rho = zeros(n,n_r);
% thick = zeros(n,n_r);
% d = zeros(1,n_r);
% Q = zeros(n,n_r);
% 
% for k = 1:n_r 
%     [y] = spsynexp(corr,n); % samples exponential sequence
%     % Needs to be changed based on relation with velocity
%     % Q having the same distribution as velocity
%     std_exp_Q = 0.3*Q_avg;
%     Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_exp_Q./std(y(:,1)));   
%     std_exp = 0.3*velavg; % Standard deviation of fluctuations 30%
%     vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_exp./std(y(:,1)));
%     while (min(vel(:,k)) <= 920 || max(vel(:,k)) >=7200 || min(Q(:,k))<=6)
%         [y] = spsynexp(corr,n);
%         vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_exp./std(y(:,1)));
%         Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_exp_Q./std(y(:,1)));
%     end
% %     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
%     std_exp_den = 0.3*den_avg;
%     rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_exp_den./std(y(:,1)));
%     thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
%     d(k) = corr.*thick_plastic; % Correlation length in terms of unit of thickness of plastic
% end

% % ---------------------- Gaussian sequence ---------------------------- 
% n_r = 100; % Number of realizations
% n = 200; % number of layers
% corr = 7; % Correlation length of Gaussian function
% vel = zeros(n,n_r);
% rho = zeros(n,n_r);
% thick = zeros(n,n_r);
% d = zeros(1,n_r);
% Q = zeros(n,n_r);
% 
% for k = 1:n_r 
%     [y] = spsyngs(corr,n); % samples Gaussian sequence
%     % Needs to be changed based on relation with velocity
%     % Q having the same distribution as velocity
%     std_gauss_Q = 0.3*Q_avg;
%     Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_gauss_Q./std(y(:,1)));   
%     std_gauss = 0.3*velavg; % Standard deviation of fluctuations 30%
%     vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_gauss./std(y(:,1)));
%     while (min(vel(:,k)) <= 920 || max(vel(:,k)) >=7200 || min(Q(:,k))<=6)
%         [y] = spsyngs(corr,n);
%         vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_gauss./std(y(:,1)));
%         Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_gauss_Q./std(y(:,1)));
%     end
% %     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
%     std_gauss_den = 0.3*den_avg;
%     rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_gauss_den./std(y(:,1)));
%     thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
%     d(k) = corr.*thick_plastic; % Correlation length in terms of unit of thickness of plastic
% end
% 


% % ---------------------- Poisson sequence ---------------------------- 
% n_r = 100; % Number of realizations
% n = 200; % number of layers
% beta = 5; % beta for Poisson's series
% d = zeros(1,n_r);
% vel = zeros(n,n_r);
% rho = zeros(n,n_r);
% thick = zeros(n,n_r);
% Q = zeros(n,n_r);
% 
% for k = 1:n_r 
%     [y,mean_thick] = spsynps(beta,n); % samples Poisson's sequence
%     std_poiss = 0.5*velavg;
%     vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_poiss./std(y(:,1)));
% %     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
%     std_poiss_den = 0.5*den_avg;
%     rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_poiss_den./std(y(:,1)));
%     thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
%     % mean thickness of Poisson media
%     d(k) = mean_thick*thick_plastic;
%     % Needs to be changed based on relation with velocity
%     % Q having the same distribution as velocity
%     std_poiss_Q = 0.5*Q_avg;
%     Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_poiss_Q./std(y(:,1)));
% end
% 

% ---------------------- Fractal sequence ---------------------------- 
n_r = 100; % Number of realizations
n = 200; % number of layers
beta = -0.8; % spectral exponent
vel = zeros(n,n_r);
rho = zeros(n,n_r);
thick = zeros(n,n_r);
d = zeros(1,n_r);
Q = zeros(n,n_r);

for k = 1:n_r 
    y = spsynfrac(beta,n); % samples using fractal sequence
    % Needs to be changed based on relation with velocity
    % Q having the same distribution as velocity
    std_frac_Q = 0.3*Q_avg;
    Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_frac_Q./std(y(:,1)));   
    std_frac = 0.3*velavg; % Standard deviation of fluctuations 30%
    vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_frac./std(y(:,1)));
    while (min(vel(:,k)) <= 920 || max(vel(:,k)) >=7200 || min(Q(:,k))<=6)
        y = spsynfrac(beta,n);
        vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_frac./std(y(:,1)));
        Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_frac_Q./std(y(:,1)));
    end
%     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
    std_frac_den = 0.3*den_avg;
    rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_frac_den./std(y(:,1)));
%     vel(:,k) = velavg; % Testing with no impedance but only Q contrast
%     rho(:,k) = den_avg; % Testing with no impedance but only Q contrast
    thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
    d(k) = thick_plastic; % Correlation length is taken as the thickness of a single layer
end


%% Calculations using Kennett-Frazer Constant Q model

lyr = zeros(n,3);
freq = [1e-3 logspace(2,9,100)];
veldisp = zeros(length(freq), n_r);
Qdisp = zeros(length(freq), n_r);
lambdaoverd = zeros(length(freq), n_r);
invq_eff_emt = zeros(length(freq), n_r);
q_rt = zeros(length(freq), n_r);


for k=1:n_r
    
    lyr(:,1) = vel(:,k);
    lyr(:,2) = rho(:,k);
    lyr(:,3) = thick(:,k);

    Q_layer = Q(:,k);

    lambdadom = 25*d(k); % Dominant wavelength
    vel_avg = mean(lyr(:,1)); % Average velocity
    fdom = vel_avg/lambdadom;


    lyr1 = repmat(lyr,100, 1);
    Q_layer1 = repmat(Q_layer, 100, 1);

    [freq,veldisp1, ~, ~,Qdisp1]=kenfdispslowQ(lyr1,freq, Q_layer1,2*pi*fdom);

    vel_emt1 = veffemt(lyr1,freq,Q_layer1,2*pi*fdom);
    
    q_emt1 = 1./(qeffemt(lyr1,freq,Q_layer1,2*pi*fdom));
    
    % Change the cutoff frequency beyond which there should be EMT based on
    % problem
%     veldisp1(freq<=10000) = vel_emt1(freq<=10000); % For Poisson medium
    veldisp1(freq<=6000) = vel_emt1(freq<=6000); % For Gaussian medium
    veldisp(:,k) = veldisp1;
%     Qdisp1(freq<=15000) = q_emt1(freq<=15000); % For Poisson medium
    Qdisp1(freq<=6000) = q_emt1(freq<=6000); % For Gaussian medium
    Qdisp(:,k) = Qdisp1;
    lambdaoverd(:,k) = (veldisp1./freq)./d(k);

    invq_eff_emt(:,k) = qeffemt(lyr1,freq,Q_layer1,2*pi*fdom);
    [~,q_rt(:,k)] = velrt_visco(lyr1,freq,Q_layer1,2*pi*fdom);
end




%% Calculations using Kennett-Frazer Constant Q model (elastic)

lyr = zeros(n,3);
freq = [1e-3 logspace(2,9,100)];
veldisp_elas = zeros(length(freq), n_r);
Qdisp_elas = zeros(length(freq), n_r);
lambdaoverd_elas = zeros(length(freq), n_r);

for k=1:n_r
    
    lyr(:,1) = vel(:,k);
    lyr(:,2) = rho(:,k);
    lyr(:,3) = thick(:,k);

    Q_layer = Q(:,k);

    lambdadom = 25*d(k); % Dominant wavelength
    vel_avg = mean(lyr(:,1)); % Average velocity
    fdom = vel_avg/lambdadom;


    lyr1 = repmat(lyr,100, 1);
    Q_layer1 = repmat(Q_layer, 100, 1);

    [freq,veldisp1_elas, ~, ~,Qdisp1_elas]=kenfdispslowQ(lyr1,freq, Q_layer1.*1e10,2*pi*fdom);

    vel_emt1_elas = veffemt(lyr1,freq,Q_layer1.*1e10,2*pi*fdom);
    
    q_emt1_elas = 1./(qeffemt(lyr1,freq,Q_layer1.*1e10,2*pi*fdom));
    
    % Change the cutoff frequency beyond which there should be EMT based on
    % problem
%     veldisp1_elas(freq<=10000) = vel_emt1_elas(freq<=10000); % For Poisson medium
    veldisp1_elas(freq<=6000) = vel_emt1_elas(freq<=6000); % For Gaussian medium
    veldisp_elas(:,k) = veldisp1_elas;
%     Qdisp1_elas(freq<=15000) = q_emt1_elas(freq<=15000); % For Poisson medium
    Qdisp1_elas(freq<=6000) = q_emt1_elas(freq<=6000); % For Gaussian medium
    Qdisp_elas(:,k) = Qdisp1_elas;
    lambdaoverd_elas(:,k) = (veldisp1_elas./freq)./d(k);
end




%% Plotting the velocity dispersion curves

figure;
semilogx(lambdaoverd,veldisp, '-', 'Color', [0.83 0.82 0.78], 'LineWidth', 3);
hold on;
semilogx(mean(lambdaoverd,2), mean(veldisp,2), '-k', 'LineWidth', 3);
mean_lambdaoverd = mean(lambdaoverd,2);
mean_veldisp = mean(veldisp,2);

semilogx(1e-1, mean_veldisp(find(mean_lambdaoverd <= 1e-1,1)), 'pk',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');

% % For Poisson medium
% semilogx(1e2, mean_veldisp(find(mean_lambdaoverd <= 1e2,1)), '^k',...
%     'MarkerSize', 15, 'MarkerFaceColor', 'k');

% For Gaussian medium
semilogx(1e3, mean_veldisp(find(mean_lambdaoverd <= 1e3,1)), '^k',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');

xlabel('\lambda/d'); 
% xlim([1e-1 1e2]); % For Poisson medium
xlim([1e-1 1e3]); % For Gaussian medium
ylabel('Velocity (m/sec)');
% ylim([1500 4000]); % For Poisson medium
ylim([2000 4500]); % For Gaussian medium



figure;
semilogx(lambdaoverd_elas,veldisp_elas, '-', 'Color', [0.83 0.82 0.78], 'LineWidth', 3);
hold on;
semilogx(mean(lambdaoverd_elas,2), mean(veldisp_elas,2), '-k', 'LineWidth', 3);
xlabel('\lambda/d'); 
% xlim([1e-1 1e2]); % For Poisson medium
xlim([1e-1 1e3]); % For Gaussian medium
ylabel('Velocity (m/sec)');
% ylim([1500 4000]); % For Poisson medium
ylim([2000 4500]); % For Gaussian medium


%% Plotting the Q curves

fig1 = figure;
set(fig1, 'Units','inches', 'Position',[0 0 20 7],'PaperPositionMode','auto');
set(gca,...
'Units','normalized',... 
'FontUnits','points',... 
'FontWeight','normal',... 
'FontSize',16,... 
'FontName','Times')


subplot 131
semilogx(lambdaoverd, 1./Qdisp, '-', 'Color', [0.83 0.82 0.78], 'LineWidth', 3);
hold on;
semilogx(mean(lambdaoverd,2), 1./harmmean(Qdisp,2), '-k', 'LineWidth', 3);
xlabel('\lambda/d'); 
% xlim([1e-1 1e2]); % For Poisson medium
xlim([1e-1 1e3]); % For Gaussian medium
ylabel('1/Q');
% ylim([0 1.5]); % For Poisson medium
ylim([0 1]); % For Gaussian medium
set(gca, 'Layer', 'top'); % To bring axis on top
title('(i) Effective (Viscoelastic)', 'FontWeight', 'bold');



% figure;
subplot 132;
semilogx(lambdaoverd_elas, 1./Qdisp_elas, '-', 'Color', [0.83 0.82 0.78], 'LineWidth', 3);
hold on;
semilogx(mean(lambdaoverd_elas,2), 1./harmmean(Qdisp_elas,2), '-k', 'LineWidth', 3);
xlabel('\lambda/d'); 
% xlim([1e-1 1e2]); % For Poisson medium
xlim([1e-1 1e3]); % For Gaussian medium
ylabel('1/Q');
% ylim([0 1.5]); % For Poisson medium
ylim([0 1]); % For Gaussian medium
set(gca, 'Layer', 'top'); % To bring axis on top
title('(ii) Eff-scattering (elastic)', 'FontWeight', 'bold');



% figure;
subplot 133;
semilogx(lambdaoverd, 1./Qdisp - 1./Qdisp_elas, '-', 'Color', [0.83 0.82 0.78], 'LineWidth', 3);
hold on;
xx = mean(lambdaoverd,2);
yy = (1./harmmean(Qdisp,2)-1./harmmean(Qdisp_elas,2));
xx_spline = logspace(log10(min(xx)),log10(max(xx)), 30);
yy_spline = spline(xx,yy,xx_spline);
semilogx(mean(lambdaoverd,2), (1./harmmean(Qdisp,2)-1./harmmean(Qdisp_elas,2)), '-k', 'LineWidth', 3);
% semilogx(xx_spline, yy_spline, '-k', 'LineWidth', 3);
xlabel('\lambda/d'); 
% xlim([1e-1 1e2]); % For Poisson medium
xlim([1e-1 1e3]); % For Gaussian medium
ylabel('1/Q');
% ylim([0 0.5]); % For Poisson medium
ylim([0 0.5]); % For Gaussian medium
set(gca, 'Layer', 'top'); % To bring axis on top
title('(iii) Effective-intrinsic', 'FontWeight', 'bold');

% % For Poisson medium
% % Finding the row index for lambdaoverd of 1e-1 and 1e2
% % Temporarily manually selected row 17 for lambdaoverd = 1e2 
% % and row 70 for lambdaoverd = 1e-1
% 
% semilogx(95, harmmean(invq_eff_emt(17,:)), '^k',...
%     'MarkerSize', 15, 'MarkerFaceColor', 'k');
% semilogx(0.15, 1./harmmean(q_rt(70,:)), 'pk',...
%     'MarkerSize', 15, 'MarkerFaceColor', 'k');
% 

% % For Gaussian medium
% % Finding the row index for lambdaoverd of 1e-1 and 1e3
% % Temporarily manually selected row 17 for lambdaoverd = 1e2 
% % and row 70 for lambdaoverd = 1e-1

% semilogx(975, harmmean(invq_eff_emt(2,:)), '^k',...
%     'MarkerSize', 15, 'MarkerFaceColor', 'k');
% semilogx(0.15, 1./harmmean(q_rt(69,:)), 'pk',...
%     'MarkerSize', 15, 'MarkerFaceColor', 'k');

% For Fractal medium
% Finding the row index for lambdaoverd of 1e-1 and 1e4
% Temporarily manually selected row 17 for lambdaoverd = 1e2 
% and row 70 for lambdaoverd = 1e-1

semilogx(975, harmmean(invq_eff_emt(2,:)), '^k',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');
semilogx(0.15, 1./harmmean(q_rt(69,:)), 'pk',...
    'MarkerSize', 15, 'MarkerFaceColor', 'k');



%% Plotting the velocity and Q models for all except Poisson medium 

realizations = [1, 50, 100];
j=0;
figure;
for i=1:length(realizations)
    subplot (3,3,j+1)
    plot(vel(:,realizations(i)),cumsum(thick(:,realizations(i))), '-k');
    axis ij; 
    xlabel('Velocity (m/sec)'); ylabel('');
    ylim([0 0.12]);

    subplot (3,3,j+2)
    plot(rho(:,realizations(i))./1e3,cumsum(thick(:,realizations(i))), '-k');
    axis ij; 
    xlabel('Density (gm/cc)'); ylabel('');
    ylim([0 0.12]);
    title(['Realization number = ' num2str(realizations(i))], 'FontWeight', 'bold');

    subplot (3,3,j+3)
%     plot(Q_layer,cumsum(lyr(:,3)), '-k');
     plot(Q(:,realizations(i)),cumsum(thick(:,realizations(i))), '-k');
    axis ij; 
    xlabel('Q'); ylabel('');
    ylim([0 0.12]);
    j=j+3;
end


%% Plot only for Poisson medium (stairstep)

figure;
subplot 131
[xx1,yy1] = stairs(vel(:,1), cumsum(thick(:,1)));
patch(xx1,yy1,'k');
axis ij; axis tight;
xlabel('Realization number 1'); ylabel('Depth (m)');
set(gca, 'XTickLabel', '');

subplot 132
[xx2,yy2] = stairs(vel(:,49), cumsum(thick(:,49)));
patch(xx2,yy2,'k');
axis ij; axis tight;
xlabel('Realization number 50'); ylabel('Depth (m)');
set(gca, 'XTickLabel', '');

subplot 133
[xx3,yy3] = stairs(vel(:,100), cumsum(thick(:,100)));
patch(xx3,yy3,'k');
axis ij; axis tight;
xlabel('Realization number 100'); ylabel('Depth (m)');
set(gca, 'XTickLabel', '');

% ---------------------------------------------------------------------------
% %% Calculation using waveforms 
% % -------------------- Using sourcewvlt function --------------------------
% % 
% lambdadom = 25*d; % Dominant wavelength
% vel_avg = mean(vel); % Average velocity
% fdom = vel_avg/lambdadom;
% % fdom =  1.6539e+04; % obtained using the lambdadom method for Gaussian
% dt_nyquist = 1/(2*fdom);
% dt = 0.0351*dt_nyquist;
% N_sample = round(dt_nyquist./dt);
% time = 0:dt:300e-3; % 10 time more time than required
% if (rem(length(time),2) ~= 0)
%     time = time(1:end-1);
% end
% wvlt = sourcewvlt;
% wvlt(end:length(time)) = 0;
% 
% % ------------------------- Waveforms ------------------------------------- 
% skip = 10; % Skipping calculation of waveforms
% skip1 = 1; % Skipping display of waveforms
% n_skip = 1:skip:n; 
% pz = zeros(length(time), n/skip);
% time_rt = zeros(1,n/skip);
% time_emt = zeros(1,n/skip);
% time_picked = zeros(1,n/skip);
% wavelength = zeros(1,n/skip);
% Q_phase = zeros(1,n/skip);
% figure;
% offset = 0;
% ytick_location = zeros(n/skip,1);
% j = 0; k = 0;
% 
% disp ('Initializing the calculations');
% 
% for i =1:skip:n
%     disp (['Iteration no ' num2str(i) '/' num2str(n)]);
%     j = j+1;
%     lyr = zeros(i,3);
%     lyr(:,1) = vel(1:i);
%     lyr(:,2) = rho(1:i);
%     lyr(:,3) = thick(1:i);
%     Q_layer = Q(1:i);
%     velinf_layer = velinf(1:i);
% 
%     % [wz,pz_temp] = kennet(lyr,wvlt,dt,2,1,-1);
%     % [wz,pz_temp] = kennettQ2(lyr,wvlt,dt,2,1,-1,Q_layer,2*pi*fdom);
%     [wz,pz_temp] = kennettQ3(lyr,wvlt,dt,2,1,-1,velinf_layer,2*pi*(fdom/2));
%     pz(:,j) = pz_temp; 
%     % Theoretical time for RT limit (elastic)
%     time_rt(j) = sum(thick(1:i)./velinf(1:i));
%     [max_index_rt] = find(time'>=time_rt(j));
%     
%     % Theoretical time for EMT limit (elastic)
%     f = thick(1:i)./sum(thick(1:i));% Fractional volumes
%     den_avg = sum(f.*rho(1:i));
%     vel_emt = sqrt(sum(f./(rho(1:i).*(vel(1:i).^2))).*den_avg);
%     vel_emt = 1./vel_emt;
%     % vel_emt = veffemt1(lyr,100,velinf_layer,2*pi*(fdom/2));
%     time_emt(j) = sum(thick(1:i))./real(vel_emt);
%     [max_index_emt] = find(time'>=time_emt(j));
%     
%     % Finding the time near 1% of the first peak
%     % [max_index] = find(time'>=time_rt(j) & pz(:,j)>=0.1*max(pz(:,j)));
%     [max_index] = find(pz(:,j)>=0.1*max(pz(:,j)));
%     time_picked(j) = time(max_index(1));
%     
%     % Finding the dominant wavelength
%     % Intepolating time and synthetic to get values closer to 0
%     time_temp = time(max_index(1)):dt/10000:time(max_index(1))+N_sample*100*dt;
%     pz_temp1 = interp1(time,pz(:,j),time_temp);
%     % Finding the zero crossings
%     [sec_index,time_zero] = crossing(pz_temp1,time_temp,0,'linear');
%     
%     velocity_cal = sum(lyr(:,3))/time(max_index(1));
%     % Dominant Wavlength
%     wavelength(i) = velocity_cal * (time_zero(2)-time(max_index(1)));
%     
%     % Plotting the synthetic seismogram
%     if (rem(i,skip1) == 0)
%         k = k+1;
%         plot((time)./1e-6,pz(:,j)+offset, 'k');
%         hold on;
%         ytick_location(k) = pz(1,j)+offset;
%         plot(time(max_index_rt(1))./1e-6, pz_temp(max_index_rt(1))+offset,'*b');
%         plot(time(max_index_emt(1))./1e-6, pz_temp(max_index_emt(1))+offset,'xb');
%         plot(time(max_index(1))./1e-6,pz_temp(max_index(1))+offset,'or');
%         plot(time_zero(2)./1e-6,pz_temp1(sec_index(2))+offset,'^c');
%     end
%     offset = offset - 1.5;
%     
%     % Spectral ratio method for Q estimation
%     % Window of 30 microsec from the first arrival time
%     % Window should be changed based on requirements
%     
%     wvlt_window = wvlt(time >=0 & time <= 30e-6);
%     pz_temp2 = pz(:,j);
%     pz_window = pz_temp2(time >= time(max_index(1)) & time <= time(max_index(1))+30e-6);
%     
%     % wvlt_window = wvlt_window';
%     
%     % taking the fourier transform of the input signals
%     [S1,f] = fftrl(wvlt_window,0:dt:(length(wvlt_window)-1)*dt); % Signal at x1
%     [S2,f] = fftrl(pz_window,time(max_index(1)):dt:time(max_index(1))+(length(pz_window)-1)*dt); % Signal at x2
%     
%     % Band limited frequency
%     % Should be changed based on source frequency
%     f = f(f>=5e4 & f<=30e4);
%     S1 = S1(f>=5e4 & f<=30e4);
%     S2 = S2(f>=5e4 & f<=30e4);
%     
%     % Taking ratio of the spectra
%     ratio_S = log(abs(S2)./abs(S1));
%     
%     % Linear least square fit through the logspectrumratio (lsr)
%     pcoeff=polyfit(f,ratio_S,1);
% 
%     Q_phase(i) = -pi*time(max_index(1))./pcoeff(1); % delta_t in this case is time of 1st arrival
%     
%     
% end
% disp ('Plotting');
% 
% xlim([0 240]);
% xlabel('Time (microseconds)'); ylabel('Number of discs');
% title('Transmitted seismogram');
% grid on;
% y1=get(gca,'ylim'); 
% set(gca,'YTick',ytick_location(end:-1:1));
% set(gca,'YTickLabel',n_skip(end:-1:1));
% 
% figure;
% plot(1:skip:n, time_rt./1e-6, '-k');
% hold on;
% plot(1:skip:n, time_emt./1e-6, '-k');
% plot(1:skip:n, time_picked./1e-6, '--k');
% xlabel('Layer number');
% ylabel('Travel time in microseconds');
% 
% figure;
% semilogx(wavelength, 1./Q_phase, 'ok');
% grid on;
% xlabel('\lambda');
% ylabel('1/Q');
% 
% 
% % Plotting the last 
% figure;
% freq = logspace(-10,10,1000);
% lyr1 = repmat(lyr,10, 1);
% Q_layer1 = repmat(Q_layer, 10, 1);
% [freq,veldisp, vel_rt, ~]=kenfdispslowQ(lyr1,freq, Q_layer1,2*pi*fdom);
% % [freq,veldisp, vel_rt, vel_emt]=kenfdispslow(lyr1,freq);
% vel_emt = veffemt(lyr1,freq,Q_layer1,2*pi*fdom);
% veldisp(freq<=1500) = vel_emt(freq<=1500);
% % veldisp(freq<=500) = vel_emt;
% semilogx((veldisp./freq)./d, veldisp, '--k');
% hold all;
% semilogx(wavelength(end)./d, velocity_cal, 'ok');
% xlim([1e-1 1e3]);
% grid on; grid minor; 
% xlabel('\lambda/d'); ylabel('Velocity'); 
% legend('Kennett-Frazer', 'From waveform at the end of layer');

