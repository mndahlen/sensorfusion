clear
close all
%% Setup
% Speed of sound
c = 343;
% Microphone positions for config 1 and 2
P1 = [-2.130, 3.546, 2.173, 3.554, 1.816, -3.458, -1.902, -3.200];
P2 = [-0.093, 3.786, -1.497, 3.765, -2.679, 2.062, -2.851, 0.176];
%% Load pulse times
load("preprocess.mat");
% Converting to distance
tphat_data = tphat_data*c;
tphat_calibration = tphat_calibration*c;
%% Calibration
% We assume the arrival time of every sensor is equal, then each sensor has
% added gaussian noise with variance sigma_k and bias b_k. That means the average is the optimal
% estimation of the arrival time which we subtract from each measurement, 
% leaving us with gaussian noise with mean 0. Each sensor is expected to 
% have its own gaussian noise variance. 

% Excluding final datapoint because its an outlier
tphat_calibration = tphat_calibration(:,1:end-1);
tphat_data = tphat_data(:,1:end-1);
N = size(tphat_calibration);
N = N(2);

% Mean arrival distance of each pulse
calibration_mean = mean(tphat_calibration,1);
% Deviation from mean arrival distance for each sensor
calibration_deviation = tphat_calibration - calibration_mean;
% First and second moment
sensor_bias = mean(calibration_deviation,2);
sensor_var = var(calibration_deviation',1);
sensor_standard_deviation = sqrt(sensor_var);

%% Setup noise distributions

% r0
r0_R1 = diag(sensor_var(1:4));
r0_R2 = diag(sensor_var(5:8));
b_1 = sensor_bias(1:4);
b_2 = sensor_bias(5:8);
r0_PE1 = ndist(b_1, r0_R1);
r0_PE2 = ndist(b_2, r0_R2);

% reference
reference_R1 = [sensor_var(2) + sensor_var(1), sensor_var(1), sensor_var(1);
               sensor_var(1), sensor_var(3) + sensor_var(1), sensor_var(1);
               sensor_var(1), sensor_var(1), sensor_var(4) + sensor_var(1)];
reference_R2 = [sensor_var(6) + sensor_var(5), sensor_var(5), sensor_var(5);
               sensor_var(5), sensor_var(7) + sensor_var(5), sensor_var(5);
               sensor_var(5), sensor_var(5), sensor_var(8) + sensor_var(5)];
b_1 = [sensor_var(2) - sensor_var(1); sensor_var(3) - sensor_var(1); sensor_var(4)- sensor_var(1)];
b_2 = [sensor_var(6) - sensor_var(5); sensor_var(7) - sensor_var(5); sensor_var(8)- sensor_var(5)];
reference_PE1 = ndist(b_1, reference_R1);
reference_PE2 = ndist(b_2, reference_R2);

%% Plot histograms
% Histfit (can put into subplot later)
% for i = 1:8
%     fig = histfit(calibration_deviation(i,:));
%     title("Noise histogram for sensor", num2str(i));
%     filename = "figs/noise_sensor_" + num2str(i) + ".png";
%     saveas(gcf, filename);
% end

%% Sensor models
S_reference_1 = sensormod('reference_tdoa', [2,0,3,8]);
S_reference_1.x0 = [0,0];
S_reference_1.pe = reference_PE1;
S_reference_1.th = P1;
S_reference_1.fs = 2;
S_r0_1 = sensormod('r0_tdoa', [3,0,4,8]);
S_r0_1.x0 = [0,0,0];
S_r0_1.pe = r0_PE1;
S_r0_1.th = P1;
S_r0_1.fs = 2;
S_reference_2 = sensormod('reference_tdoa', [2,0,3,8]);
S_reference_2.x0 = [0,0];
S_reference_2.pe = reference_PE2;
S_reference_2.th = P2;
S_reference_2.fs = 2;
S_r0_2 = sensormod('r0_tdoa', [3,0,4,8]);
S_r0_2.x0 = [0,0,0];
S_r0_2.pe = r0_PE2;
S_r0_2.th = P2;
S_r0_2.fs = 2;

%% Configuration analysis
if 0
    resolution = 0.1;
    xlow=-4;
    xhigh=4;
    
    % Config 1
    figure()
    y = simulate(S_reference_1,0);
    plot(S_reference_1)
    hold on;
    lh2(S_reference_1,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh))
    axis([xlow,xhigh,xlow,xhigh])
    crlb(S_reference_1,y)
    view([-90,90])
    
    % Config 2
    figure()
    y = simulate(S_reference_2,0);
    plot(S_reference_2)
    hold on;
    lh2(S_reference_2,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh))
    axis([xlow,xhigh,xlow,xhigh])
    crlb(S_reference_2,y)
    view([-90,90])
end
% Config 1 is better

%% Localization
% I will do a and d for configuration 1
% Remove sensor bias
tphat_data = tphat_data - sensor_bias;
% config 1
tphat_data = tphat_data(1:4,:);

% a)
if 0       
    P_inv = inv(r0_R1);
    num_points = 100;
    resolution_xy = (max(P1)-min(P1))/num_points;

    num_points_r = 10000;
    low_r = 0;
    high_r = 100;
    eps_r = 200; %0.5*343 = 171.5
    eps_r_2 = -140;

    estimates = zeros(3,132);
    for time = [1:132]
        time
        resolution_r = (high_r-low_r)/num_points;
        min_loss = 9999999999999;
        min_x = 0;
        min_y = 0;
        min_r = 0;
        for x = [min(P1):resolution_xy:max(P1)]
            for y = [min(P1):resolution_xy:max(P1)]
                for r = [low_r:resolution_r:high_r]
                    y_hat = r0_tdoa(0,[x;y;r],0,P1);
                    y_real = tphat_data(:,time);
                    loss = (y_real-y_hat)'*P_inv*(y_real-y_hat);
                    if loss < min_loss
                        min_x = x;
                        min_y = y;
                        min_r = r;
                        min_loss = loss;
                    end
                end
            end
        end
        low_r = min_r-eps_r_2;
        high_r = min_r+eps_r;
        estimates(:,time) = [min_x;min_y;min_r];
    end
    save estimates_a estimates
end
    
% d)
if 0
    y_permute = [-1 1 0 0; 
                 -1 0 1 0; 
                 -1 0 0 1];
    tphat_data_p = y_permute*tphat_data;
    P_inv = inv(reference_R1);
    num_points = 1000;
    resolution_xy = (max(P1)-min(P1))/num_points;

    estimates = zeros(2,132);
    for time = [1:132]
        time
        min_loss = 9999999999999;
        min_x = 0;
        min_y = 0;
        for x = [min(P1):resolution_xy:max(P1)]
            for y = [min(P1):resolution_xy:max(P1)]
                y_hat = reference_tdoa(0,[x;y],0,P1);
                y_real = tphat_data_p(:,time);
                loss = (y_real-y_hat)'*P_inv*(y_real-y_hat);
                if loss < min_loss
                    min_x = x;
                    min_y = y;
                    min_loss = loss;
                end
            end
        end
        [min_x, min_y, min_loss]
        estimates(:,time) = [min_x;min_y];
    end
    save estimates_d estimates
end

%% Tracking
load estimates_d.mat
cv2d = exmotion('cv2d');
ctcv2d = exmotion('ctcv2d');

% a) Location sensor model

% Motion model 1
if 0
    z = sig(estimates');
    z.fs=2;
    S_location_1 = sensormod(@(t,x,u,th)[x(1,:); x(2,:)], [2,0,2,0]);
    S_location_1.pe = trace(reference_PE1.P)/3*eye(2);
    S_location_1.fs = 2;
    S_location_1_cv2d = addsensor(cv2d, S_location_1);
    x = ekf(S_location_1_cv2d, z);
    xplot2(x, 'conf', 90);
    view([-90 90])
    saveas(gcf, "figs/tracking_ekf_motion_1","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_1","png");
end

% Motion model 2
if 0
    S_location_1 = sensormod(@(t,x,u,th)[x(1,:); x(2,:)], [2,0,2,0]);
    S_location_1.pe = trace(reference_PE1.P)/3*eye(2);
    S_location_1.fs = 2;
    S_location_1_ctcv2d = addsensor(ctcv2d, S_location_1);
    x = ekf(S_location_1_ctcv2d, z);
    xplot2(x, 'conf', 90);
    view([-90 90])
    saveas(gcf, "figs/tracking_ekf_motion_2","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_2","png");
end

% b) Reference TDOA sensor model
y_reference_tdoa = zeros(3,132);
for time =[1:132]
    y_reference_tdoa(:,time) = reference_tdoa(0,[estimates(1,time),estimates(2,time)],0,P1);
end

% Motion model 1
if 1
    z = sig(y_reference_tdoa');
    z.fs=2;
    z.Py = reference_PE1.P;
    S_reference_1_cv2d = addsensor(cv2d, S_reference_1);
    S_reference_1_cv2d.x0 = [1,1,1,1]';
    [x,V] = ekf(S_reference_1_cv2d, z)
end
%% 7 Uncertainty
% How much can I perturb before really bad (maximum dist)
