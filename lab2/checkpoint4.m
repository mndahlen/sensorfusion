% Remove gyroscope?
% Mag MU floating around?
% Outlier rejection: Use mean? (Turn off in phone)
% Outlier rejection: Do nothing (?) Do not perform the MU.

if 1
    load saved/acc_mean.mat
    load saved/acc_cov.mat
    load saved/gyr_mean.mat
    load saved/gyr_cov.mat
    load saved/mag_mean.mat
    load saved/mag_cov.mat
    
    calAcc.m = acc_mean;
    calAcc.R = acc_cov;
    calGyr.m = gyr_mean;
    calGyr.R = gyr_cov;
    calMag.m = mag_mean;
    calMag.R = mag_cov;
    
    [xhat, meas] = checkpoint4_EKF('', calAcc, calGyr, calMag)
    
    save saved/xhat_checkpoint_4_side xhat 
    save saved/meas_checkpoint_4_side meas
end

% FLAT
if 0
    % Estimates vs Google baseline - Flat start
    load saved/xhat_checkpoint_4_flat.mat
    load saved/meas_checkpoint_4_flat.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), meas.orient(1, idxs), color="blue", DisplayName="x 1");
    plot(meas.t(:,idxs), meas.orient(2, idxs), color="red",  DisplayName="x 2");
    plot(meas.t(:,idxs), meas.orient(3, idxs), color="black",  DisplayName="x 3");
    plot(meas.t(:,idxs), meas.orient(4, idxs), color="green",  DisplayName="x 4");
    plot(meas.t(:, idxs), xhat.x(1, idxs), color="blue", LineStyle="--", DisplayName="xh 1");
    plot(meas.t(:, idxs), xhat.x(2, idxs), color="red", LineStyle="--", DisplayName="xh 2");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="black", LineStyle="--", DisplayName="xh 3");
    plot(meas.t(:, idxs), xhat.x(4, idxs), color="green", LineStyle="--", DisplayName="xh 4");
    legend
    title("Estimates vs Google - Flat")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_flat_estimate.png");
    hold off
    
    % Error Estimates vs Google baseline - Flat start
    load saved/xhat_checkpoint_4_flat.mat
    load saved/meas_checkpoint_4_flat.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), xhat.x(1, idxs) - meas.orient(1, idxs), color="blue", DisplayName="dx 1");
    plot(meas.t(:,idxs), xhat.x(2, idxs) - meas.orient(2, idxs), color="red",  DisplayName="dx 2");
    plot(meas.t(:,idxs), xhat.x(3, idxs) - meas.orient(3, idxs), color="black",  DisplayName="dx 3");
    plot(meas.t(:,idxs), xhat.x(4, idxs) - meas.orient(4, idxs), color="green",  DisplayName="dx 4");
    legend
    title("Error - Flat")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_flat_error.png");
    hold off
end

% SIDE
if 0
    % Estimates vs Google baseline - Side start
    load saved/xhat_checkpoint_4_side.mat
    load saved/meas_checkpoint_4_side.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), meas.orient(1, idxs), color="blue", DisplayName="x 1");
    plot(meas.t(:,idxs), meas.orient(2, idxs), color="red",  DisplayName="x 2");
    plot(meas.t(:,idxs), meas.orient(3, idxs), color="black",  DisplayName="x 3");
    plot(meas.t(:,idxs), meas.orient(4, idxs), color="green",  DisplayName="x 4");
    plot(meas.t(:, idxs), xhat.x(1, idxs), color="blue", LineStyle="--", DisplayName="xh 1");
    plot(meas.t(:, idxs), xhat.x(2, idxs), color="red", LineStyle="--", DisplayName="xh 2");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="black", LineStyle="--", DisplayName="xh 3");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="green", LineStyle="--", DisplayName="xh 4");
    legend
    title("Estimates vs Google - Side")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_side_estimate.png");
    hold off
    
    % Error Estimates vs Google baseline - Side start
    load saved/xhat_checkpoint_4_side.mat
    load saved/meas_checkpoint_4_side.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), xhat.x(1, idxs) - meas.orient(1, idxs), color="blue", DisplayName="dx 1");
    plot(meas.t(:,idxs), xhat.x(2, idxs) - meas.orient(2, idxs), color="red",  DisplayName="dx 2");
    plot(meas.t(:,idxs), xhat.x(3, idxs) - meas.orient(3, idxs), color="black",  DisplayName="dx 3");
    plot(meas.t(:,idxs), xhat.x(4, idxs) - meas.orient(4, idxs), color="green",  DisplayName="dx 4");
    legend
    title("Error - Side")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_side_error.png");
    hold off
end

% FLAT outlier
if 0
    % Estimates vs Google baseline - Flat start
    load saved/xhat_checkpoint_4_flat_outlier.mat
    load saved/meas_checkpoint_4_flat_outlier.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), meas.orient(1, idxs), color="blue", DisplayName="x 1");
    plot(meas.t(:,idxs), meas.orient(2, idxs), color="red",  DisplayName="x 2");
    plot(meas.t(:,idxs), meas.orient(3, idxs), color="black",  DisplayName="x 3");
    plot(meas.t(:,idxs), meas.orient(4, idxs), color="green",  DisplayName="x 4");
    plot(meas.t(:, idxs), xhat.x(1, idxs), color="blue", LineStyle="--", DisplayName="xh 1");
    plot(meas.t(:, idxs), xhat.x(2, idxs), color="red", LineStyle="--", DisplayName="xh 2");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="black", LineStyle="--", DisplayName="xh 3");
    plot(meas.t(:, idxs), xhat.x(4, idxs), color="green", LineStyle="--", DisplayName="xh 4");
    legend
    title("Estimates vs Google - Flat")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_flat_estimate_outlier.png");
    hold off
    
    % Error Estimates vs Google baseline - Flat start
    load saved/xhat_checkpoint_4_flat_outlier.mat
    load saved/meas_checkpoint_4_flat_outlier.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), xhat.x(1, idxs) - meas.orient(1, idxs), color="blue", DisplayName="dx 1");
    plot(meas.t(:,idxs), xhat.x(2, idxs) - meas.orient(2, idxs), color="red",  DisplayName="dx 2");
    plot(meas.t(:,idxs), xhat.x(3, idxs) - meas.orient(3, idxs), color="black",  DisplayName="dx 3");
    plot(meas.t(:,idxs), xhat.x(4, idxs) - meas.orient(4, idxs), color="green",  DisplayName="dx 4");
    legend
    title("Error - Flat")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_flat_error_outlier.png");
    hold off
end

% SIDE outlier
if 0
    % Estimates vs Google baseline - Side start
    load saved/xhat_checkpoint_4_side_outlier.mat
    load saved/meas_checkpoint_4_side_outlier.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), meas.orient(1, idxs), color="blue", DisplayName="x 1");
    plot(meas.t(:,idxs), meas.orient(2, idxs), color="red",  DisplayName="x 2");
    plot(meas.t(:,idxs), meas.orient(3, idxs), color="black",  DisplayName="x 3");
    plot(meas.t(:,idxs), meas.orient(4, idxs), color="green",  DisplayName="x 4");
    plot(meas.t(:, idxs), xhat.x(1, idxs), color="blue", LineStyle="--", DisplayName="xh 1");
    plot(meas.t(:, idxs), xhat.x(2, idxs), color="red", LineStyle="--", DisplayName="xh 2");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="black", LineStyle="--", DisplayName="xh 3");
    plot(meas.t(:, idxs), xhat.x(3, idxs), color="green", LineStyle="--", DisplayName="xh 4");
    legend
    title("Estimates vs Google - Side")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_side_estimate_outlier.png");
    hold off
    
    % Error Estimates vs Google baseline - Side start
    load saved/xhat_checkpoint_4_side_outlier.mat
    load saved/meas_checkpoint_4_side_outlier.mat
    meas_size = size(meas.t)
    idxs = 2:meas_size(2) 
    fig = figure();
    hold on
    plot(meas.t(:,idxs), xhat.x(1, idxs) - meas.orient(1, idxs), color="blue", DisplayName="dx 1");
    plot(meas.t(:,idxs), xhat.x(2, idxs) - meas.orient(2, idxs), color="red",  DisplayName="dx 2");
    plot(meas.t(:,idxs), xhat.x(3, idxs) - meas.orient(3, idxs), color="black",  DisplayName="dx 3");
    plot(meas.t(:,idxs), xhat.x(4, idxs) - meas.orient(4, idxs), color="green",  DisplayName="dx 4");
    legend
    title("Error - Side")
    xlabel("t [s]")
    ylabel("values")
    saveas(gcf, "figs/checkpoint_4_side_error_outlier.png");
    hold off
end
