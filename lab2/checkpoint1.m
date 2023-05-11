% [xhat, meas] = filterTemplate()
% 
% save meas_checkpoint_1 meas

load meas_checkpoint_1.mat

%% Histogram of all measurements (acc, gyr, mag)
if 0
    % acc 1
    fig = histfit(meas.acc(1, ~any(isnan(meas.acc), 1)));
    title("acc axis 1");
    filename = "figs/acc_1.png";
    saveas(gcf, filename);
    
    % acc 2
    fig = histfit(meas.acc(2, ~any(isnan(meas.acc), 1)));
    title("acc axis 2");
    filename = "figs/acc_2.png";
    saveas(gcf, filename);
    
    % acc 3
    fig = histfit(meas.acc(3, ~any(isnan(meas.acc), 1)));
    title("acc axis 3");
    filename = "figs/acc_3.png";
    saveas(gcf, filename);
    
    % gyr 1
    fig = histfit(meas.gyr(1, ~any(isnan(meas.gyr), 1)));
    title("gyr axis 1");
    filename = "figs/gyr_1.png";
    saveas(gcf, filename);
    
    % gyr 2
    fig = histfit(meas.gyr(2, ~any(isnan(meas.gyr), 1)));
    title("gyr axis 2");
    filename = "figs/gyr_2.png";
    saveas(gcf, filename);
    
    % gyr 3
    fig = histfit(meas.gyr(3, ~any(isnan(meas.gyr), 1)));
    title("gyr axis 3");
    filename = "figs/gyr_3.png";
    saveas(gcf, filename);
    
    % mag 1
    fig = histfit(meas.mag(1, ~any(isnan(meas.mag), 1)));
    title("mag axis 1");
    filename = "figs/mag_1.png";
    saveas(gcf, filename);
    
    % mag 2
    fig = histfit(meas.mag(2, ~any(isnan(meas.mag), 1)));
    title("mag axis 2");
    filename = "figs/mag_2.png";
    saveas(gcf, filename);
    
    % mag 3
    fig = histfit(meas.mag(3, ~any(isnan(meas.mag), 1)));
    title("mag axis 3");
    filename = "figs/mag_3.png";
    saveas(gcf, filename);
end

%% Mean and variance of all measurements
if 1
    % acc
    acc_mean = mean(meas.acc(:, ~any(isnan(meas.acc), 1)), 2)
    acc_cov = cov((meas.acc(:, ~any(isnan(meas.acc), 1)))')
    
    % gyr
    gyr_mean = mean(meas.gyr(:, ~any(isnan(meas.gyr), 1)), 2)
    gyr_cov = cov((meas.gyr(:, ~any(isnan(meas.gyr), 1)))')
    
    % mag
    mag_mean = mean(meas.mag(:, ~any(isnan(meas.mag), 1)), 2)
    mag_cov = cov((meas.mag(:, ~any(isnan(meas.mag), 1)))')

    save acc_mean acc_mean
    save acc_cov acc_cov
    save gyr_mean gyr_mean
    save gyr_cov gyr_cov
    save mag_mean mag_mean
    save mag_cov mag_cov
end

%% Trends of all measurements
if 0
    fig = plot(meas.t(:,~any(isnan(meas.acc), 1)), meas.acc(:, ~any(isnan(meas.acc), 1)));
    title("Acceleration");
    filename = "figs/acc_plot.png";
    saveas(gcf, filename);
    
    fig = plot(meas.t(:,~any(isnan(meas.gyr), 1)), meas.gyr(:, ~any(isnan(meas.gyr), 1)));
    title("Gyroscope");
    filename = "figs/gyr_plot.png";
    saveas(gcf, filename);
    
    fig = plot(meas.t(:,~any(isnan(meas.mag), 1)), meas.mag(:, ~any(isnan(meas.mag), 1)));
    title("Magnetometer");
    filename = "figs/mag_plot.png";
    saveas(gcf, filename);
end