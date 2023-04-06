%% Load recorded data (both audio recording and positional recordings).
% Make sure to add the data directory to the path
filename_data = 'data';
filename_calibration = 'calibration';
[data, fs] = audioread(sprintf("data/%s.wav", filename_data));
[calibration, fs] = audioread(sprintf("data/%s.wav", filename_calibration));
load(sprintf("data/%s.mat", filename_data));
load(sprintf("data/%s.mat", filename_calibration));
%% Setup
% Channels to detect pulses in
channels = 1:8;
% Number of pulses per second
npt = 2;
% Plot during processing
plot_on = 0;
% Type of pulse - chirp or ofdm
type = 'chirp';
% Minimal prominence of detected peaks, tuneable
min_prominence = 0.0;
%% Detect pulse times and indices
tphat_data = SFlabFindPulseTimes(data, channels, fs, ...
                                    npt, plot_on, type, min_prominence);
tphat_calibration = SFlabFindPulseTimes(calibration, channels, fs, ...
                                    npt, plot_on, type, min_prominence);

%% Calibration
% We assume the arrival time of every sensor is equal, then each sensor has
% added gaussian noise with mean 0. That means the average is the optimal
% estimation of the arrival time which we subtract from each measurement, 
% leaving us with gaussian noise with mean 0. Each sensor is expected to 
% have its own gaussian noise variance, which 

% Excluding final datapoint because its an outlier
tphat_calibration = tphat_calibration(:,1:end-1);
N = size(tphat_calibration);
N = N(2);

% Mean arrival time of each pulse
calibration_mean = mean(tphat_calibration,1);

% Deviation from mean arrival time for each sensor
calibration_deviation = tphat_calibration - calibration_mean;

% Deviations give measurement variance of each sensor 
% (unbiased estimate of variance)
calibration_variance = sum(calibration_deviation.^2,2)/(N-1);
calibration_standard_deviation = sqrt(calibration_variance);

% Histfit (can put into subplot later)
f1 = histfit(calibration_deviation(1,:));
f2 = histfit(calibration_deviation(2,:));
f3 = histfit(calibration_deviation(3,:));
f4 = histfit(calibration_deviation(4,:));
f5 = histfit(calibration_deviation(5,:));
f6 = histfit(calibration_deviation(6,:));
f7 = histfit(calibration_deviation(7,:));
f8 = histfit(calibration_deviation(8,:));
%% Sensor models
% https://www.youtube.com/watch?v=h2hZK8LTRSs

% Model 1
h1 = sensormod

%% Microphone locations (find in Qualisys -- in m)
mic_locations = [-2.130,  3.546;
                  2.173,  3.554;
                  1.816, -3.458;
                 -1.902, -3.200;
                 -0.093,  3.786;
                 -1.497,  3.765;
                 -2.679,  2.062;
                 -2.851,  0.176]';
x0 = position(1, 2:end);
position = position'; % Transpose position as ROS saves as N x M