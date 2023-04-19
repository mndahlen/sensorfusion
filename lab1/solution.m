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

% Reference
R_1 = diag(sensor_var(1:4));
R_2 = diag(sensor_var(5:8));
b_1 = sensor_bias(1:4);
b_2 = sensor_bias(5:8);
reference_PE1 = ndist(b_1, R_1);
reference_PE2 = ndist(b_2, R_2);

% Residual
R_1 = [sensor_var(2) + sensor_var(1), sensor_var(1), sensor_var(1);
       sensor_var(1), sensor_var(3) + sensor_var(1), sensor_var(1);
       sensor_var(1), sensor_var(1), sensor_var(4) + sensor_var(1)];
R_2 = [sensor_var(6) + sensor_var(5), sensor_var(5), sensor_var(5);
       sensor_var(5), sensor_var(7) + sensor_var(5), sensor_var(5);
       sensor_var(5), sensor_var(5), sensor_var(8) + sensor_var(5)];
b_1 = [sensor_var(2) - sensor_var(1); sensor_var(3) - sensor_var(1); sensor_var(4)- sensor_var(1)];
b_2 = [sensor_var(6) - sensor_var(5); sensor_var(7) - sensor_var(5); sensor_var(8)- sensor_var(5)];
residual_PE1 = ndist(b_1, R_1);
residual_PE2 = ndist(b_2, R_2);

%% Plot histograms
% Histfit (can put into subplot later)
% for i = 1:8
%     fig = histfit(calibration_deviation(i,:));
%     title("Noise histogram for sensor", num2str(i));
%     filename = "figs/noise_sensor_" + num2str(i) + ".png";
%     saveas(gcf, filename);
% end

%% Sensor models
S1_1 = sensormod('residual_tdoa', [2,0,3,8]);
S1_1.pe = residual_PE1;
S1_1.th = P1;
S1_1.fs = 2;
S2_1 = sensormod('reference_tdoa', [3,0,4,8]);
S2_1.pe = reference_PE1;
S2_1.th = P1;
S2_1.fs = 2;
S1_2 = sensormod('residual_tdoa', [2,0,3,8]);
S1_2.pe = residual_PE2;
S1_2.th = P2;
S1_2.fs = 2;
S2_2 = sensormod('reference_tdoa', [3,0,4,8]);
S2_2.pe = reference_PE2;
S2_2.th = P2;
S2_2.fs = 2;

%% Configuration analysis
if 0
    resolution = 0.1;
    xlow=-4;
    xhigh=4;
    
    % Config 1
    figure()
    y = simulate(S1_1,0);
    plot(S1_1)
    hold on;
    lh2(S1_1,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh))
    axis([xlow,xhigh,xlow,xhigh])
    crlb(S1_1,y)
    view([-90,90])
    
    % Config 2
    figure()
    y = simulate(S1_2,0);
    plot(S1_2)
    hold on;
    lh2(S1_2,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh))
    axis([xlow,xhigh,xlow,xhigh])
    crlb(S1_2,y)
    view([-90,90])
end
% Config 1 is better

%% Localization
if 0
    % I will do a and d
    
    % Remove sensor bias
    tphat_data = tphat_data - sensor_bias;
    % config 1
    tphat_data = tphat_data(1:4,:)
    
    % a)
    sigobj = sig(tphat_data);
    sigobj.fs = 2;
    
    
    % b)
    y_permute = [-1 1 0 0; 
                 -1 0 1 0; 
                 -1 0 0 1];
    tphat_data_p = y_permute*tphat_data
    sigobj = sig(tphat_data_l);
    sigobj.fs = 2;
end