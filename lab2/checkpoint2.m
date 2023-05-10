
load acc_mean acc_mean
load acc_cov acc_cov
load gyr_mean gyr_mean
load gyr_cov gyr_cov
load mag_mean mag_mean
save mag_cov mag_cov

[xhat, meas] = checkpoint2_EKF()