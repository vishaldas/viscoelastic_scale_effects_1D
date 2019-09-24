close all
clear all
% Creating the samples
vel_plastic = 2487; % 2487m/sec
vel_steel = 5535; % 5535m/sec
den_plastic = 1.210*1000; % 1.210 g/cc
den_steel = 7.900 *1000; % 7.900 g/cc
D = 52e-3; % Length of sample 52e-3
P1 = 0.5; % Proportion of plastic
P2 = 1-P1; % Proportion of steel
Q_plastic =10; % 10
Q_steel = 20; % 20

fdom = 200e3;
dt_nyquist = 1/(2*fdom);
dt = 0.0351*dt_nyquist;
N_sample = round(dt_nyquist./dt);
time = 0:dt:70e-3-dt;

fdom1 = 115e3;

% -------------------- Periodicity ----------------------------------------
M_period = [1 2 3 5 6 7 8 9 10 12 14 16 32 64 128 256]; % Periodicity
% M_period = [1 2 3 4 5 6 7 8 9 10 12 14 16 32 64];
% M_period = [1 4 10 16 32 128]; % Periodicity

velocity_cal = zeros(length(M_period),1);
wavelength = zeros(length(M_period),1);
R = zeros(length(M_period),1);
Q_phase = zeros(length(M_period),1);
tf = cell(length(M_period),1);
kf = cell(length(M_period),1);
vel_eff = cell(length(M_period),1);
invq_eff = cell(length(M_period),1);

pz_pick = zeros(length(M_period),1);
time_pick = zeros(length(M_period),1);


offset = 0;
ytick_location = zeros(length(M_period),1);

disp ('Initializing the calculations');

% For elastic media

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
    
    
    % Transfer functions for range of frequency
    om = 2.*pi.*logspace(3,6,1000);
    
%     [wz,pz_temp,tf_temp] = kennet(lyr,wvlt,dt,2,1,-1);
    [tf_temp] = kennettQ2_tf(lyr,om,2,0,Q_layer*1e10,2*pi*fdom1);

    tf(i) = {tf_temp};
        
end

disp ('Plotting');


% For viscoelastic media

tf_ve = cell(length(M_period),1);

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
    
    
    % Transfer functions for range of frequency
    om = 2.*pi.*logspace(3,6,1000);
    
%     [wz,pz_temp,tf_temp] = kennet(lyr,wvlt,dt,2,1,-1);
    [tf_temp] = kennettQ2_tf(lyr,om,2,0,Q_layer,2*pi*fdom1);

    tf_ve(i) = {tf_temp};
        
end

disp ('Plotting');

%% Plotting few of the power spectrum and comparing elastic and viscoelastic
% Q total will be related to the width of the peak at half power
% High Q = narrow peak

figure;
plot_M = [1 2 7 12 13 16];

for i = 1:length(plot_M)
    tf1 = cell2mat(tf(plot_M(i)));
    tf2 = cell2mat(tf_ve(plot_M(i)));
    power_tf1 = (abs(tf1(2:end,3))).^2;
    power_rf1 = (abs(tf1(2:end,2))).^2;
    
    power_tf2 = (abs(tf2(2:end,3))).^2;
    power_rf2 = (abs(tf2(2:end,2))).^2;
        
    tf1_freq = tf1(2:end,1);
    tf2_freq = tf2(2:end,1);
    
    subplot (3,2,i);
    semilogx(tf1_freq,power_tf1, '-k');
    hold on;
    semilogx(tf2_freq,power_tf2, '-r');
    semilogx(tf1_freq, power_rf1, '--k');
    semilogx(tf2_freq, power_rf2, '--r');
    set(gca, 'Layer', 'top');
    set(gca, 'XTick', logspace(3,6,4));
%     y1=get(gca,'ylim');
%     semilogx([fdom fdom], y1, '--k');
    xlim([1e3 1e5]);
    ylim([0 1]);
    xlabel('Frequency'); ylabel('Power');
    title(['Number of periods = ' num2str(M_period(plot_M(i)))]); 
end
legend('Transmission', 'Reflection', 'Transmission', 'Reflection');


% %% Plotting power spectrum (subplot) 
% 
% figure;  
%   
% 
% for i = 1:length(M_period)
%     tf1 = cell2mat(tf(i));
%     power_tf1 = (abs(tf1(2:end,3))).^2;
%     power_rf1 = (abs(tf1(2:end,2))).^2;
%     
%     power_total = power_tf1 + power_rf1;
%     
%     % Calculate phase from transfer functions 
%     phase_tf = unwrap(angle(tf1(2:end,3)));
%     
%     tf1_freq = tf1(2:end,1);
%     
%     subplot (4,4,i);
%     ax = gca;
%     semilogx(tf1_freq,power_tf1, '-k');
%     ylabel('Power');
%     hold on;
% %     yyaxis right;
%     semilogx(tf1_freq, power_rf1, '-b');
% %     ax.YScale = 'log';
% %     ax.YColor = 'b';
%     ylabel('Power');
%     ylim([0 1]);
% 
% %     y1=get(gca,'ylim');
% %     semilogx([fdom fdom], y1, '--k');
%     xlabel('Frequency'); 
%     title(['Number of periods = ' num2str(M_period(i))]); 
% end
% legend('Transmission', 'Reflection');
% 
% 
% %% Plotting (full scale) power spectrum
% 
% 
% for i = 1:length(M_period)
%     tf1 = cell2mat(tf(i));
%     power_tf1 = (abs(tf1(2:end,3))).^2;
%     power_rf1 = (abs(tf1(2:end,2))).^2;
%         
%     tf1_freq = tf1(2:end,1);
%     figure;
%     loglog(tf1_freq,power_tf1, '-k');
%     hold on;
%     loglog(tf1_freq, power_rf1, '-b');
%     set(gcf,'pos',[10 10 812 650]);
% %     y1=get(gca,'ylim');
% %     semilogx([fdom fdom], y1, '--k');
%     xlabel('Frequency'); ylabel('power');
%     title(['Number of periods = ' num2str(M_period(i))]);
%     legend('Transmission', 'Reflection', 'Position', [0.716338263185349 0.926410258427644 0.176108370640595 0.0746153825979965]);
%     set(gca, 'Layer', 'top'); % To bring axis on top
% end
% 
% %% Plotting few of the power spectrum 
% % Q total will be related to the width of the peak at half power
% % High Q = narrow peak
% 
% figure;
% plot_M = [1 2 7 12 13 16];
% 
% for i = 1:length(plot_M)
%     tf1 = cell2mat(tf(plot_M(i)));
%     power_tf1 = (abs(tf1(2:end,3))).^2;
%     power_rf1 = (abs(tf1(2:end,2))).^2;
%         
%     tf1_freq = tf1(2:end,1);
%     subplot (3,2,i);
%     semilogx(tf1_freq,power_tf1, '-k');
%     hold on;
%     semilogx(tf1_freq, power_rf1, '--k');
%     set(gca, 'Layer', 'top');
%     set(gca, 'XTick', logspace(3,6,4));
% %     y1=get(gca,'ylim');
% %     semilogx([fdom fdom], y1, '--k');
%     xlim([1e3 1e5]);
%     ylim([0 1]);
%     xlabel('Frequency'); ylabel('Power');
%     title(['Number of periods = ' num2str(M_period(plot_M(i)))]); 
% end
% legend('Transmission', 'Reflection');
% 
% 
% 
% 
% %% Approx lambda/d for each of the selected plots
% 
% lambdaoverd = [0.3426 0.6922 2.0451 2.1868 5.3811 44.7974];
% 
% M_period = [1 2 8 16 32 256]; % Periodicity
% 
% 
% title(['Number of periods = ' num2str(M_period(i)) ', \lambda/d  \approx ' ...
%     num2str(lambdaoverd(i))]);
% 
% %% Relative phase shift with respect to fastest (ray theory) 
% % and then the phase lag for different periods
% 
% % phase_tf = zeros(length(om)-1, length(M_period));
% % phase_rf = zeros(length(om)-1, length(M_period));
% phase_tf1 = zeros(length(om)-1, length(M_period));
% phase_rf1 = zeros(length(om)-1, length(M_period));
% 
% for i = 1:length(M_period)
%     tf1 = cell2mat(tf(i));
%     tf_relative = cell2mat(tf(1));
%     
% %     % Calculate phase from transfer functions 
% %     phase_tf(:,i) = (angle(tf1(2:end,3)));
% %     phase_rf(:,i) = (angle(tf1(2:end,2)));
% 
%     % Relative phase difference
%     
%     phase_tf1(:,i) = angle(tf1(2:end,3) ./tf_relative(2:end,3));
%     phase_rf1(:,i) = angle(tf1(2:end,2) ./tf_relative(2:end,2));
% 
%     
%     tf1_freq = tf1(2:end,1);
%     
% end
% 
% % % Relative phase shift
% % phase_lag_tf = phase_tf(:, 1:end) - ...
% %     repmat(phase_tf(:,1), [1 length(M_period)]);
% % phase_lag_rf = phase_rf(:, 1:end) - ...
% %     repmat(phase_rf(:,1), [1 length(M_period)]);
% 
% figure;
% 
% for i = 1:length(M_period)
%     subplot (4,4,i);
%     semilogx(tf1_freq,phase_tf1(:,i), '-k');
%     hold on;
%     semilogx(tf1_freq,phase_rf1(:,i), '-b');
%     xlabel('Frequency'); ylabel('Phase shift (in radians)');
%     title(['Number of periods = ' num2str(M_period(i))]); 
%     axis tight;
%     ylim([-pi pi]);
%     
% end
% legend('Transmission', 'Reflection');
% 
% 
% figure;
% plot_M = [1 2 7 12 13 16];
% for i = 1:length(plot_M)
%     subplot (3,2,i);
%     semilogx(tf1_freq,phase_tf1(:,plot_M(i)), '-k');
%     hold on;
%     semilogx(tf1_freq,phase_rf1(:,plot_M(i)), '--k');
%     xlabel('Frequency'); ylabel('Phase shift (in radians)');
%     title(['Number of periods = ' num2str(M_period(plot_M(i)))]); 
%     axis tight;
%     ylim([-pi pi]);
%     xlim([1e3 1e5]);
%     
% end
% legend('Transmission', 'Reflection');
