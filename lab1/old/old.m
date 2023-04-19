%% Setup
% Speed of sound
c = 343;
% Microphone positions for config 1 and 2
P = [-2.130,  3.546;
       2.173,  3.554;
       1.816, -3.458;
      -1.902, -3.200;
      -0.093,  3.786;
      -1.497,  3.765;
      -2.679,  2.062;
      -2.851,  0.176];
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
calibration_deviation_bias = mean(calibration_deviation,2)
% Deviations give measurement variance of each sensor 
% (unbiased estimate of variance)
calibration_variance = sum(calibration_deviation.^2,2)/(N-1);
calibration_standard_deviation = sqrt(calibration_variance);

% Histfit (can put into subplot later)
for i = 1:8
    fig = histfit(calibration_deviation(i,:));
    title("Noise histogram for sensor", num2str(i));
    filename = "figs/noise_sensor_" + num2str(i) + ".png";
    saveas(gcf, filename);
end

R_1 = diag(calibration_variance(1:4));
R_2 = diag(calibration_variance(5:8));

%% Sensor models
% The .m file
run("init_sensor_models.m")

%% Config analysis NLS
% Since the car drove inside the convex hull of the microphones
% we can use it as the bounds for grid-search
min_x = min(P(:,1)) - 0.1;
max_x = max(P(:,1)) + 0.1;
min_y = min(P(:,2)) - 0.1;
max_y = max(P(:,2)) + 0.1;
precision = 0.01;
size_x = round((max_x-min_x)/precision);
size_y = round((max_y-min_y)/precision);
NLS_1 = zeros(size_y, size_x);
NLS_2 = zeros(size_y, size_x);

% Setup residuals and variance
y= tphat_data(:,1);
res_1 = [y(2)-y(1),y(3)-y(1),y(4)-y(1),y(2)-y(3),y(4)-y(3),y(2)-y(4)]';
res_2 = [y(6)-y(5),y(7)-y(5),y(8)-y(5),y(6)-y(7),y(8)-y(7),y(6)-y(8)]';
R_res_1 = diag([R_1(2,2)+R_1(1,1),R_1(3,3)+R_1(1,1),R_1(4,4)+R_1(1,1),R_1(2,2)+R_1(3,3),R_1(4,4)+R_1(3,3),R_1(2,2)+R_1(4,4)]);
R_res_1_inv = inv(R_res_1)/100;
R_res_2 = diag([R_2(2,2)+R_2(1,1),R_2(3,3)+R_2(1,1),R_2(4,4)+R_2(1,1),R_2(2,2)+R_2(3,3),R_2(4,4)+R_2(3,3),R_2(2,2)+R_2(4,4)]);
R_res_2_inv = inv(R_res_2)/100;

T = size_x*size_y;
i = 1;
% Precompute for faster eval
A1 = res_1'*R_res_1_inv*res_1;
B1 = res_1'*R_res_1_inv;
C1 = R_res_1_inv*res_1;
A2 = res_2'*R_res_2_inv*res_2;
B2 = res_2'*R_res_2_inv;
C2 = R_res_2_inv*res_2;

for idx_x = 1:size_x
    for idx_y = 1:size_y
        if mod(i,10000)==0
            display(i/T)
        end
        i = i + 1;
        x = [idx_x*precision + min_x, idx_y*precision + min_y];
        % Calc h for config 1
        h_1 = [norm(x-P(2,:)) - norm(x-P(1,:)), 
               norm(x-P(3,:)) - norm(x-P(1,:)), 
               norm(x-P(4,:)) - norm(x-P(1,:)), 
               norm(x-P(2,:)) - norm(x-P(3,:)), 
               norm(x-P(4,:)) - norm(x-P(3,:)), 
               norm(x-P(2,:)) - norm(x-P(4,:))];
        % Calc h for config 2
        h_2 = [norm(x-P(6,:)) - norm(x-P(5,:)), 
               norm(x-P(7,:)) - norm(x-P(5,:)), 
               norm(x-P(8,:)) - norm(x-P(5,:)), 
               norm(x-P(6,:)) - norm(x-P(7,:)), 
               norm(x-P(8,:)) - norm(x-P(7,:)), 
               norm(x-P(6,:)) - norm(x-P(8,:))];
        NLS_1(idx_y,idx_x) = A1 - B1*h_1 - h_1'*C1 + h_1'*R_res_1_inv*h_1;
        NLS_2(idx_y,idx_x) = A2 - B2*h_2 - h_2'*C2 + h_2'*R_res_2_inv*h_2;
    end
end
% Plot
% Config 1xxx
[value, index] = min(NLS_1(:));
[row, col] = ind2sub(size(NLS_1), index);
imagesc(NLS_1)
colorbar;
hold on;
plot(col,row,'r*');
text(col, row, num2str(min_x + col*precision) + "," + num2str(min_y + row*precision));
for i = [1,2,3,4]
    plot((P(i,1)-min_x)/precision,(P(i,2)-min_y)/precision,'bx');
end
title("Config 1 NLS loss at k=1");
filename = "figs/config_1_loss_analysis.png";
saveas(gcf, filename);
hold off;

% Config 2
[value, index] = min(NLS_2(:));
[row, col] = ind2sub(size(NLS_2), index);
imagesc(NLS_2)
colorbar; 
hold on;
plot(col,row,'r*');
text(col, row, num2str(min_x + col*precision) + "," + num2str(min_y + row*precision));
for i = [5,6,7,8]
    plot((P(i,1)-min_x)/precision,(P(i,2)-min_y)/precision,'bx');
end
title("Config 2 NLS loss at k=1");
filename = "figs/config_2_loss_analysis.png";
saveas(gcf, filename);
hold off;