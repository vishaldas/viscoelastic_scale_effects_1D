%% Generating the random medium

close all;
clear all;

addpath('/Volumes/GoogleDrive/My Drive/Research/Codes/srbtools');

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
% 

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


% ---------------------- Poisson sequence ----------------------------
n_r = 100; % Number of realizations
n = 200; % number of layers
beta = 5; % beta for Poisson's series
d = zeros(1,n_r);
vel = zeros(n,n_r);
rho = zeros(n,n_r);
thick = zeros(n,n_r);
Q = zeros(n,n_r);

for k = 1:n_r
    [y,mean_thick] = spsynps(beta,n); % samples Poisson's sequence
    std_poiss = 0.5*velavg;
    vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_poiss./std(y(:,1)));
%     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
    std_poiss_den = 0.5*den_avg;
    rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_poiss_den./std(y(:,1)));
    thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
    % mean thickness of Poisson media
    d(k) = mean_thick*thick_plastic;
    % Needs to be changed based on relation with velocity
    % Q having the same distribution as velocity
    std_poiss_Q = 0.5*Q_avg;
    Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_poiss_Q./std(y(:,1)));
end


% % ---------------------- Fractal sequence ----------------------------
% n_r = 100; % Number of realizations
% n = 200; % number of layers
% beta = -0.8; % spectral exponent
% vel = zeros(n,n_r);
% rho = zeros(n,n_r);
% thick = zeros(n,n_r);
% d = zeros(1,n_r);
% Q = zeros(n,n_r);
% 
% for k = 1:n_r
%     y = spsynfrac(beta,n); % samples using fractal sequence
%     % Needs to be changed based on relation with velocity
%     % Q having the same distribution as velocity
%     std_frac_Q = 0.3*Q_avg;
%     Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_frac_Q./std(y(:,1)));
%     std_frac = 0.3*velavg; % Standard deviation of fluctuations 30%
%     vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_frac./std(y(:,1)));
%     while (min(vel(:,k)) <= 920 || max(vel(:,k)) >=7200 || min(Q(:,k))<=6)
%         y = spsynfrac(beta,n);
%         vel(:,k) = velavg + (y(:,1)-mean(y(:,1))).*(std_frac./std(y(:,1)));
%         Q(:,k) = Q_avg + (y(:,1)-mean(y(:,1))).*(std_frac_Q./std(y(:,1)));
%     end
%     %     rho(:,k) = (1.741.*(vel(:,k)./1000).^0.25).*1000; % Gardener's relation
%     std_frac_den = 0.3*den_avg;
%     rho(:,k) = den_avg + (y(:,1)-mean(y(:,1))).*(std_frac_den./std(y(:,1)));
%     %     vel(:,k) = velavg; % Testing with no impedance but only Q contrast
%     %     rho(:,k) = den_avg; % Testing with no impedance but only Q contrast
%     thick(:,k) = linspace(thick_plastic, thick_plastic, length(vel(:,k))).'; % Constant thickness
%     d(k) = thick_plastic; % Correlation length is taken as the thickness of a single layer
% end

%% Calculating the transfer functions

lyr = zeros(n,3);
tf = cell(n_r,1);
tf_ve = cell(n_r,1);

% Transfer functions for range of frequency
om = 2.*pi.*logspace(3,6,1000);


for k = 1:n_r
    
    lyr(:,1) = vel(:,k);
    lyr(:,2) = rho(:,k);
    lyr(:,3) = thick(:,k);
    
    Q_layer = Q(:,k);
    
    lambdadom = 25*d(k); % Dominant wavelength
    vel_avg = mean(lyr(:,1)); % Average velocity
    fdom = vel_avg/lambdadom;
    
    lyr1 = repmat(lyr,100, 1);
    Q_layer1 = repmat(Q_layer, 100, 1);
    
    
    %     [wz,pz_temp,tf_temp] = kennet(lyr,wvlt,dt,2,1,-1);
    [tf_temp] = kennettQ2_tf(lyr,om,2,0,Q_layer,2*pi*fdom);
    [tf_temp1] = kennettQ2_tf(lyr,om,2,0,Q_layer.*1e10,2*pi*fdom);
    
    tf_ve(k) = {tf_temp};
    tf(k) = {tf_temp1};
end


%% Plotting

close all;
figure;

iter = [1 50 100];

for k = 1:length(iter)
    tf1 = cell2mat(tf(iter(k)));
    tf2 = cell2mat(tf_ve(iter(k)));
    power_tf1 = (abs(tf1(2:end,3))).^2;
    power_rf1 = (abs(tf1(2:end,2))).^2;
    
    power_tf2 = (abs(tf2(2:end,3))).^2;
    power_rf2 = (abs(tf2(2:end,2))).^2;
    
    power_total = power_tf1 + power_rf1;
    
    tf1_freq = tf1(2:end,1);
    tf2_freq = tf2(2:end,1);
    
    subplot(3,1,k);
    semilogx(tf1_freq,power_tf1, '-k');
    hold on;
    semilogx(tf2_freq,power_tf2, '-r');
    semilogx(tf1_freq, power_rf1, '--k');
    semilogx(tf2_freq, power_rf2, '--r');
    set(gca, 'Layer', 'top');
    set(gca, 'XTick', logspace(3,6,4));
    %     xlim([1e3 1e5]);
    ylim([0 1]);
    xlabel('Frequency'); ylabel('Power');
    title(['Realization number = ' num2str(iter(k))]);
    
end

legend('Trans Elas', 'Trans VE', 'Refl Elas', 'Refl VE');
