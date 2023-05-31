load DATAFILE
% checkpoint 1
load saved/meas_checkpoint_1.mat
xhat2 = xhat;
meas2 = meas;
save DATAFILE -append xhat2 meas2

% checkpoint 2
load saved/xhat_checkpoint_2_flat.mat
load saved/meas_checkpoint_2_flat.mat
xhat3 = xhat;
meas3 = meas;
save DATAFILE -append xhat3 meas3

% checkpoint 3
load saved/xhat_checkpoint_3_flat.mat
load saved/meas_checkpoint_3_flat.mat
xhat4 = xhat;
meas4 = meas;
save DATAFILE -append xhat4 meas4

% checkpoint 3 outlier
load saved/xhat_checkpoint_3_flat_outlier.mat
load saved/meas_checkpoint_3_flat_outlier.mat
xhat5 = xhat;
meas5 = meas;
save DATAFILE -append xhat5 meas5

% checkpoint 4
load saved/xhat_checkpoint_4_flat.mat
load saved/meas_checkpoint_4_flat.mat
xhat6 = xhat;
meas6 = meas;
save DATAFILE -append xhat6 meas6

% checkpoint 5
load saved/xhat_checkpoint_4_flat_outlier.mat
load saved/meas_checkpoint_4_flat_outlier.mat
xhat7 = xhat;
meas7 = meas;
save DATAFILE -append xhat7 meas7

% checkpoint 6
load saved/xhat_checkpoint_5_side.mat
load saved/meas_checkpoint_5_side.mat
xhat8 = xhat; 
meas8 = meas;
save DATAFILE -append xhat8 meas8