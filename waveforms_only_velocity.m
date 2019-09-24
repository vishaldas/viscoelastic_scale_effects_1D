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
Q_plastic =10;
Q_steel = 20;

% % ----------------- Using Ricker SRBTools as source -----------------------
%
% fdom = 200e3; % Dominant frequency
% [wvlt,dt_wvlt] = ricker(fdom);
%
% wvlt = wvlt./max(wvlt);
%
% N_sample = 10;
% dt = dt_wvlt/N_sample;
%
% time = 0:dt:70e-3;
%
% wvlt(end:length(time)) = 0;
%
% time_shift = time(wvlt == max(wvlt)); % Converting to zero phase

% % ----------------- Using sine wave as source -----------------------
%
% fdom = 200e3; % Dominant frequency
% dt_nyquist = 1/(2*fdom);
% N_sample = 100;
% dt = dt_nyquist/N_sample;
%
% time = 0:dt:70e-3;
%
% wvlt = sin(2.*pi.*fdom.*time(1:2*N_sample));
%
% wvlt(end:length(time)) = 0;

% -------------------- Using sourcewvlt function --------------------------

fdom = 200e3;
dt_nyquist = 1/(2*fdom);
dt = 0.0351*dt_nyquist;
N_sample = round(dt_nyquist./dt);
time = 0:dt:70e-3-dt;

wvlt = sourcewvlt;
wvlt(end:length(time)) = 0;

fdom1 = 115e3;

% -------------------- Periodicity ----------------------------------------
% M_period = [1 2 3 5 6 7 8 9 10 12 14 16 32 64 128 256]; % Periodicity
M_period = [1 2 3 4 5 6 7 8 9 10 12 14 16 32 64];
% M_period = [1 4 10 16 32 128]; % Periodicity

pz = zeros(length(M_period),length(wvlt));
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
    
    % Theoretical ray theory and effective medium theory velocities
    
    % Calculated using sample with two layers
    
    if (M == 1)
        
        thickness1 = lyr(1,3);
        total_thickness = sum(lyr(:,3));
        
        density1 = lyr(1,2);
        density2 = lyr(2,2);
        
        velocity_1 = lyr(1,1);
        velocity_2 = lyr(2,1);
        
        fraction_1 = thickness1/total_thickness;
        fraction_2 = 1-fraction_1;
        
        den_avg = (fraction_1*density1) + (fraction_2*density2);
        
        vel_rt = ((fraction_1/velocity_1)+(fraction_2/velocity_2)).^-1;
        
        temp1 = fraction_1/(density1*velocity_1^2) + fraction_2/ ...
            (density2*velocity_2^2);
        vel_emt = sqrt(1/(temp1*den_avg));
        
        % Theoretical ray theory and emt first arrival times
        time_rt = total_thickness/vel_rt;
        time_emt = total_thickness/vel_emt;
        
    end
    
    % Synthetic seismogram
    
%     [wz,pz_temp,tf_temp] = kennet(lyr,wvlt,dt,2,1,-1);
    [wz,pz_temp,tf_temp] = kennettQ2(lyr,wvlt,dt,2,1,-1,Q_layer,2*pi*fdom1);

    tf(i) = {tf_temp};
    pz(i,:) = pz_temp;
    % pz(i,:) = pz_temp;
    
    % Calculating K-F viscoelastic curves
    lyr_kf = repmat(lyr,10,1);
    Q_layer_kf = repmat(Q_layer, 10,1);
    if (log10(min(tf_temp(:,1))) <-6)
        freq_kf = logspace(-6,log10(max(tf_temp(:,1))), 100);
    else
        freq_kf = logspace(log10(min(tf_temp(:,1))),log10(max(tf_temp(:,1))), 100);
    end
    [freq_kf,vel_disp,~,~,Q_disp] = kenfdispslowQ(lyr_kf,freq_kf,Q_layer_kf,2*pi*fdom);
    kf(i) = {[freq_kf.', vel_disp.', Q_disp.']};
    [freq_kf,vel_disp] = kenfdispslow(lyr_kf,freq_kf);
    kf(i) = {[freq_kf, vel_disp]};
    
    % Calculating the EMT curves
    vel_eff(i) = {veffemt(lyr,freq_kf,Q_layer,2*pi*fdom)};
    invq_eff(i) = {qeffemt(lyr,freq_kf,Q_layer,2*pi*fdom)};
    
    
    % Plotting the synthetic seismogram
    
%     plot((time)./1e-6,pz(i,:)+offset, '-k');
%     wigplot_edit(pz(i,1:1000).'+offset,time(1)./1e-6,(time(2)-time(1))./1e-6,0,1.5);
    
    ytick_location(i) = pz(i,1)+offset;
    
    % Finding the time near 5% of the maximum amplitude
    [~,max_index] = find(pz(i,:)>=0.01*max(pz(i,:)));
    % [~,max_index] = find(time>=time_rt & pz(i,:)>=0.01*max(pz(i,:)));
    
    % [pk,time_peak] = findpeaks(pz(j,:),'NPeaks',1,'SortStr','descend');
    % [pk,time_peak] = findpeaks(pz(j,:),'MinPeakHeight',0.01*max(pz(j,:)));
    
%     hold on;
%     plot(time(max_index(1))./1e-6,pz_temp(max_index(1))+offset,'ok', 'MarkerSize', 10);
    
    time_pick(i) = time(max_index(1));
    pz_pick(i) = pz_temp(max_index(1));    
    
    % Calculating the velocity
    velocity_cal(i) = sum(lyr(:,3))/time(max_index(1));
    
%     % Time from first arrival pick to second zero crossing
%     
    % Intepolating time and synthetic to get values closer to 0
    time_temp = time(max_index(1)):dt/10000:time(max_index(1))+N_sample*100*dt;
    pz_temp1 = interp1(time,pz(i,:),time_temp);
    time_temp_wvlt = 0:dt/10000:N_sample*100*dt;
    wvlt_temp1 = interp1(time,wvlt,time_temp_wvlt);
%     
%     % pz_temp1 = pz_temp1./max(pz_temp1);
%     
%     % Finding the zero crossings
    [sec_index,time_zero] = crossing(pz_temp1,time_temp,0,'linear');
    [sec_index_wvlt,time_zero_wvlt] = crossing(wvlt_temp1,time_temp_wvlt,-0.01,'linear');
%     
%     % Wavlength
%     wavelength(i) = velocity_cal(i) * (time_zero(2)-time(max_index(1)));
%     %   Checking the zero crossing
%     plot(time_zero(2)./1e-6,pz_temp1(sec_index(2))+offset,'*k');
    
    % Wavelength calculation based on source frequency
    wavelength(i) = velocity_cal(i)./fdom;
    R(i) = wavelength(i)./d;
    
        
    %  -------------- Spectral ratio method taking the fmax ---------------
%     % Windowing the signal from first arrival to the second zero crossing
%     wvlt_window = wvlt(time>=0 & time<=time_zero_wvlt(4));
%     pz_temp2 = pz(i,:);
%     pz_window = pz_temp2(time>=time(max_index(1)) & time <=time_zero(3));
%     % pz_window1 = [zeros(1,max_index(1)) pz_window];
%     
%     wvlt_window = wvlt_window';
%     % Calculating the amplitude spectra
%     % taking the fourier transform of the input signals
%     
%     %      pz_window = repmat(pz_window,[1 10]);
%     %      wvlt_window = repmat(wvlt_window,[1 10]);
%     
%     
%     if (length(pz_window)>= length(wvlt_window))
%         wvlt_window(end:length(pz_window)) =0;
%     else
%         pz_window(end:length(wvlt_window)) =0;
%     end
%     
%     time_window = 0:dt:(length(wvlt_window)-1)*dt;
%     
%     [S1,f] = fftrl(wvlt_window,time_window); % Signal at x1
%     [S2,f] = fftrl(pz_window,time_window); % Signal at x2
%     
%         % Calculations only for dominant frequency
%         fmax = f(abs(S1)==max(abs(S1)));
%         S1 = S1(f==fmax);
%         S2 = S2(f==fmax);
%     % Band limited frequency
% %     % Should be changed based on source frequency
% %     f = f(f>=5 & f<=250e3);
% %     S1 = S1(f>=5 & f<=250e3);
% %     S2 = S2(f>=5 & f<=250e3);
%     
%     
%     x1 = 0; % top of the sample
%     x2 = sum(lyr(:,3)); % sample size
%     
% %     % Calculating the phase velocity
% %     % phase_S1 = unwrap(atan2(imag(S1),real(S1)));
% %     phase_S1 = unwrap(angle(S1));
% %     phase_S2 = unwrap(angle(S2));
% %     
% %     phase_diff = phase_S2-phase_S1;
% %     %
% %     %     v_phase = abs((2.*pi.*f.*(x2-x1))./(phase_diff));
%     
%     
%     % Taking ratio of the spectra
%     ratio_S = log(abs(S2)./abs(S1));
%     
%     % To be used when calculation only made for fmax
%     Q_phase(i) = (-pi.*fmax.*(x2-x1))./(velocity_cal(i).*ratio_S);
%     % To be used when calculation is done on all the frequencies
% %     Q_phase_temp = (-pi.*f.*(x2-x1))./(velocity_cal(i).*ratio_S);
% %     Q_phase(i) = nanmean(Q_phase_temp);
    
    
    % ----------------- Spectral ratio method using lsq fitting -----------
%         % Constant window method
%         % Window of 70 microsec from the first arrival time
%         % Window should be changed based on requirements
%     
%         wvlt_window = wvlt(time >=0 & time <= 70e-6);
%         pz_temp2 = pz(i,:);
%         pz_window = pz_temp2(time >= time(max_index(1)) & time <= time(max_index(1))+70e-6);
%     
%         wvlt_window = wvlt_window';
%         % taking the fourier transform of the input signals
%         [S1,f] = fftrl(wvlt_window,0:dt:(length(wvlt_window)-1)*dt); % Signal at x1
%         [S2,f] = fftrl(pz_window,time(max_index(1)):dt:time(max_index(1))+(length(pz_window)-1)*dt); % Signal at x2
%     
    
    
%         % Windowing the signal from first arrival to the second zero crossing
%         wvlt_window = wvlt(time>=0 & time<=time_zero_wvlt(4));
%         pz_temp2 = pz(i,:);
%         pz_window = pz_temp2(time>=time(max_index(1)) & time <=time_zero(3));
%     
%         wvlt_window = wvlt_window';
%     
%         if (length(pz_window)>= length(wvlt_window))
%             wvlt_window(end:length(pz_window)) =0;
%         else
%             pz_window(end:length(wvlt_window)) =0;
%         end
%         
%         wvlt_window(end:end+1000) = 0;
%         pz_window(end:end+1000) = 0;
%     
%         time_window = 0:dt:(length(wvlt_window)-1)*dt;
%         
%         % Testing the windowing of the wavelet and the seismogram
% %         subplot 122
% %         plot(time_window./1e-6, pz_window+offset, '-r');
% %         hold on
% %         plot(time_window./1e-6, wvlt_window+offset, '-b');
%     
%         % taking the fourier transform of the input signals
%         [S1,f] = fftrl(wvlt_window,time_window); % Signal at x1
%         [S2,f] = fftrl(pz_window,time_window); % Signal at x2
%     
%         % Band limited frequency
%         % 20% on both sides of the dominant frequency
%         fband = f(f>= fdom.*(1-0.2) & f<=fdom.*(1+0.2));
%         S1band = S1(f>= fdom.*(1-0.2) & f<=fdom.*(1+0.2));
%         S2band = S2(f>= fdom.*(1-0.2) & f<=fdom.*(1+0.2));
%     
%         % Taking ratio of the spectra
%         ratio_S = log(abs(S2band)./abs(S1band));
%     
%         % Linear least square fit through the logspectrumratio (lsr)
%         % pcoeff=polyfit(fband,ratio_S,1);
%         options=optimset('algorithm','active-set');
%         pcoeff = lsqlin(fband',ratio_S',[],[],0,0,[],[], [],options);
%     
%         % Method1 
% %         Q_phase(i) = -pi*time(max_index(1))./pcoeff(1); % delta_t in this case is time of 1st arrival
%         
% %         % Method2
%         x1 = 0; % top of the sample
%         x2 = sum(lyr(:,3)); % sample size
%         Q_phase(i) = -pi.*(x2-x1)./(velocity_cal(i).* pcoeff);
    
    offset = offset - 0.5;
    
end

%%
disp ('Plotting');

figure;
wigplot_edit(-pz(:,1:1000).',time(1)./1e-6,(time(2)-time(1))./1e-6,0,0.5, pz_pick.', time_pick.'./1e-6);
set(gca, 'ydir', 'reverse');
xlim([0 70]);
xlabel('Time (microseconds)'); 
% ylabel('Periodicity');
title('Transmitted seismograms');
grid on;

% Plotting the theoretical limits

y1=get(gca,'ylim');
hold on;
plot([time_rt/1e-6 time_rt/1e-6],y1, '--k');
plot([time_emt/1e-6 time_emt/1e-6],y1, '--k');
% 
% % Setting the labels
% 
% set(gca,'YTick',ytick_location(end:-1:1));
% set(gca,'YTickLabel',M_period(end:-1:1));
set(gca,'YTickLabel','');


%% This is to plot the velocity dispersion curves 
load('working_workspace.mat');

figure;
semilogx(R, velocity_cal, 'ok');
hold on;
semilogx(lambdaoverd_kf, velcalc_kf, '-*k');
semilogx(lambdaoverd_emt, vel_eff_emt, '--k');
xlabel('\lambda/d'); ylabel('Velocity (m/sec)');
legend('Velocity from waveforms','K-F velocities', 'Viscoelastic EMT');
grid on; box on; 
xlim([1e-1 1e2]); ylim([1000 5000]);


%%

% % Calculating the velocity dispersion curve
% figure;
% semilogx(R, velocity_cal, 'ok');
% grid on;
% xlabel('\lambda/d');
% ylabel('Velocity (m/sec)');
% x1lim = get(gca,'xlim');
% hold on;
% plot(x1lim,[vel_rt vel_rt], '--k');
% plot(x1lim,[vel_emt vel_emt], '--k');
% ylim([vel_emt-100 vel_rt+100]);


% figure;
% semilogx(R,(1./Q_phase), 'ok');
% grid on;
% xlabel('\lambda/d');
% ylabel('1/Q');
% ylim([0 max(1./Q_phase)+0.2]);

% %% Plotting K-F viscoelastic curves and the phase velocity for each of the samples
% 
% figure;
% 
% for i = 1:length(M_period)-1
%     kf1 = cell2mat(kf(i));
%     tf1 = cell2mat(tf(i));  
%     vel_eff_emt = cell2mat(vel_eff(i)).';
%     
%     % Calculate phase velocities from transfer functions 
%     phase_tf = unwrap(angle(tf1(:,3)));
%     vel_tf = 2.*pi.*tf1(:,1).*D./phase_tf;
%     
%     % Cleaning data below 15% of the dominant frequency
%      tf1_freq = tf1(:,1);
%      vel_tf(tf1_freq<=0.15*fdom) = nan;
%      
%      % Cleaning k-f data for emt values
%      kf1_vel = kf1(:,2);
%      kf1_freq = kf1(:,1);
%      kf1_vel(kf1_freq<=0.15*fdom) = vel_eff_emt(kf1_freq<=0.15*fdom);
%     
%     subplot (4,4,i);
%     semilogx(tf1_freq,vel_tf, '-');
%     hold all;
%     semilogx(kf1_freq, kf1_vel, '--');
%     semilogx(kf1(:,1), vel_eff_emt, '-.');
%     y11lim = get(gca,'ylim');
%     semilogx([200e3 200e3], [0 6000], '--k');
%     xlim([1e3 1e7]);
%     ylim([0 6000]);
%     xlabel('Frequency'); ylabel('Velocity');
%     title(['Periodicity = ' num2str(M_period(i))]);
%     grid on; box on; 
%     
% end
% 
%     legend('From transfer function', 'From Kennett-Frazer viscoelastic', 'EMT viscoelastic', 'Dominant Freq = 200KHz');


%% Velocity dispersion 

% figure;
% 
% freqcalc_tf = zeros(1,length(M_period));
% velcalc_tf = zeros(1,length(M_period));
% freqcalc_kf = zeros(1,length(M_period));
% velcalc_kf = zeros(1,length(M_period));
% qcalc_kf = zeros(1,length(M_period));
% 
% dcalc = D./M_period;
% 
% for i = 1:length(M_period)
%     kf1 = cell2mat(kf(i));
%     tf1 = cell2mat(tf(i));
%     vel_eff_emt = cell2mat(vel_eff(i)).';
%     
%     % Calculate phase velocities from transfer functions 
%     phase_tf = unwrap(angle(tf1(:,3)));
%     vel_tf = 2.*pi.*tf1(:,1).*D./phase_tf;
%     
%     % Cleaning data below 15% of the dominant frequency
%      tf1_freq = tf1(:,1);
%      vel_tf(tf1_freq<=0.15*fdom) = nan;
%      
%      % Cleaning k-f data for emt values
%      kf1_vel = kf1(:,2);
%      kf1_freq = kf1(:,1);
%      kf1_vel(kf1_freq<=2e3) = vel_eff_emt(kf1_freq<=2e3);
%      kf1_q = kf1(:,3);
%      
%      % Finding the velocity value at fdom 
%      temp1 = abs(tf1_freq-fdom);
%      [indx indx] = min(temp1);
%      freq_temp = logspace(log10(tf1_freq(indx-3)),log10(tf1_freq(indx+3)),10);
%      vel_temp = interp1(tf1_freq, vel_tf, freq_temp);
%      [indx1 indx1] = min(abs(freq_temp-fdom));
%      
%      freqcalc_tf(i) = freq_temp(indx1);
%      velcalc_tf(i) = vel_temp(indx1);
%      
%      temp2 = abs(kf1_freq-fdom);
%      [indx2 indx2] = min(temp2);
%      freq_temp = logspace(log10(kf1_freq(indx2-2)),log10(kf1_freq(indx2+2)),100);
%      vel_temp = interp1(kf1_freq, kf1_vel, freq_temp);
%      [indx3 indx3] = min(abs(freq_temp-fdom));
%      q_temp = interp1(kf1_freq, kf1_q, freq_temp);
%      [indx4 indx4] = min(abs(freq_temp-fdom));
%      
%      freqcalc_kf(i) = freq_temp(indx3);
%      velcalc_kf(i) = vel_temp(indx3);
%      qcalc_kf(i) = q_temp(indx4);
%      
% end

%%

% % Values of K_F corrected for stacks with small periodicity
% velcalc_kf(M_period <=2) = velcalc_tf(M_period <=2);
%      
% subplot 121
% semilogx(freqcalc_tf, velcalc_tf, 'ok');
% hold on;
% semilogx(freqcalc_kf, velcalc_kf, 'ob');
% xlim([1e5 1e6]);
% xlabel('Frequency'); ylabel('Velocity');
% legend('Transfer function of transmissivity','Kennett-Frazer viscoelastic');
% 
% lambdaoverd_tf = (velcalc_tf./freqcalc_tf)./dcalc;
% lambdaoverd_kf = (velcalc_kf./freqcalc_kf)./dcalc;
% emt_input = cell2mat(kf(17));
% kf1_freq_emt = emt_input(:,1);
% vel_eff_emt = cell2mat(vel_eff(17)).';
% invq_eff_emt = cell2mat(invq_eff(17)).';
% lambdaoverd_emt = (vel_eff_emt./kf1_freq_emt)./dcalc(17);
% lambdaoverd_emt(lambdaoverd_emt <= 5) = nan;
% subplot 122
% semilogx(lambdaoverd_tf, velcalc_tf, '--b');
% hold on;
% semilogx(R, velocity_cal, 'ok');
% semilogx(lambdaoverd_emt, vel_eff_emt, '--k');
% semilogx(lambdaoverd_kf, velcalc_kf, 'or');
% xlabel('\lambda/d'); ylabel('Velocity');
% legend('K-F velocity dispersion','Waveforms', 'Viscoelastic EMT');
% xlim([1e-1 1e2]);



%% Q values plot
% load('Q_total.mat');
% 
% figure;
% semilogx(R1,(1./Q_phase1),'xk', 'MarkerSize',10);
% hold all;
% semilogx(R,(1./Q_phase), 'ok', 'MarkerSize', 10);
% semilogx(R1, 1./Q_phase1 - 1./Q_phase, 'dk', 'MarkerSize', 10);
% semilogx(lambdaoverd_emt1, invq_eff_emt1, '--k');
% grid on;
% xlabel('\lambda/d');
% ylabel('1/Q');
% ylim([0 0.8]);
% xlim([1e-1 1e2]); 
% legend('Viscoelastic','Elastic','Intrinsic', 'Effective medium theory');
% 


%% Plotting the values used in the phase velocity calculation for 1 case

% i =15;
% kf1 = cell2mat(kf(i)).';
% tf1 = cell2mat(tf(i));
% vel_eff_emt = cell2mat(vel_eff(i)).';
% 
% % Calculate phase velocities from transfer functions
% phase_tf = unwrap(angle(tf1(:,3)));
% vel_tf = 2.*pi.*tf1(:,1).*D./phase_tf;
% 
% figure;
% subplot 231
% semilogx(tf1(:,1), vel_tf);
% xlabel('Frequency'); ylabel('Velocity');
% 
% subplot 232
% semilogx(tf1(:,1), phase_tf);
% xlabel('Frequency'); ylabel('Phase');
% 
% subplot 233
% semilogx(tf1(:,1), real(tf1(:,3)));
% hold all;
% semilogx(tf1(:,1), imag(tf1(:,3)));
% xlabel('Frequency'); ylabel('Transmissivity Transfer function');
% legend('Real part', 'Imaginary part');
% 
% tf1_freq = tf1(:,1);
% vel_tf(tf1_freq<=0.15*fdom) = nan;
% phase_tf(tf1_freq<=0.15*fdom) = nan;
% 
% subplot 234
% semilogx(tf1_freq, vel_tf);
% xlabel('Frequency'); ylabel('Velocity');
% 
% subplot 235
% semilogx(tf1_freq, phase_tf);
% xlabel('Frequency'); ylabel('Phase');


%% Q estimation from Kennett-Frazer 

