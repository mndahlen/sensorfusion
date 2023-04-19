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
%% Adjust
% Excluding final datapoint because its an outlier
tphat_calibration = tphat_calibration(:,1:end-1);

%% Save
save preprocess.mat tphat_data tphat_calibration
