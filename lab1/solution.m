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
S_reference_1 = sensormod(@reference_tdoa, [2,0,3,8]);
S_reference_1.x0 = [0,0];
S_reference_1.pe = reference_PE1;
S_reference_1.th = P1;
S_reference_1.fs = 2;
S_r0_1 = sensormod('r0_tdoa', [3,0,4,8]);
S_r0_1.x0 = [0,0,0];
S_r0_1.pe = r0_PE1;
S_r0_1.th = P1;
S_r0_1.fs = 2;
S_reference_2 = sensormod(@reference_tdoa, [2,0,3,8]);
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

    y = simulate(S_reference_1,0);
    % Plot config 1 setup and CRLB
    if 0
        figure()
%         lh2(S_reference_1,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh));
        plot(S_reference_1)
        hold on;
        axis([xlow,xhigh,xlow,xhigh])
        crlb(S_reference_1,y)
        view([-90,90]);
        ax = gca;
        ax.FontSize = 16; 
%         title("Configuration 1 and CRLB with target at [0,0]")
        saveas(gcf, "figs/config_1_setup_analysis","epsc");
        saveas(gcf, "figs/config_1_setup_analysis","png");
        hold off;
    end

    % Plot config 1 loss
    if 0
        [lh_1,x1,x2,px,px0,X1,X2] = lh2(S_reference_1,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh));
        [value, index] = min(lh_1(:));
        [row, col] = ind2sub(size(lh_1), index);
        imagesc(lh_1)
        hold on
    %     colorbar;
        plot(col,row,'r*');
        text(col, row, num2str(xlow + col*resolution) + "," + num2str(xlow + row*resolution));
        for i = [1,2,3,4]
            plot((P1(2*(i-1)+1)-xlow)/resolution,(P1(2*(i-1)+2)-xlow)/resolution,'bx');
        end
        view([90,-90]);
        ax = gca;
        ax.FontSize = 16; 
        xlabel("x1")
        ylabel("x2")
        xticks([1 20 40 60 80 ])
        xticklabels({xlow + 0, xlow+20*resolution, xlow+40*resolution, xlow+60*resolution, xlow+80*resolution})
        yticks([1 20 40 60 80 ])
        yticklabels({xlow + 0, xlow+20*resolution, xlow+40*resolution, xlow+60*resolution, xlow+80*resolution})
%         title("Configuration 1 and square loss with target at [0,0]")
        saveas(gcf, "figs/config_1_loss_analysis","epsc");
        saveas(gcf, "figs/config_1_loss_analysis","png");
        hold off
    end

    % Config 2
    y = simulate(S_reference_2,0);
    % Plot config 2 setup and CRLB
    if 0
        figure()
%         lh2(S_reference_2,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh));
        plot(S_reference_2)
        hold on
        axis([xlow,xhigh,xlow,xhigh])
        crlb(S_reference_2,y)
        view([-90,90]);
        ax = gca;
        ax.FontSize = 16; 
%         title("Configuration 2 and CRLB with target at [0,0]")
        saveas(gcf, "figs/config_2_setup_analysis","epsc");
        saveas(gcf, "figs/config_2_setup_analysis","png");
        hold off
    end

    % Plot config 2 loss
    if 0
        [lh_2,x1,x2,px,px0,X1,X2] = lh2(S_reference_2,y,(xlow:resolution:xhigh),(xlow:resolution:xhigh));
        [value, index] = min(lh_2(:));
        [row, col] = ind2sub(size(lh_2), index);
        imagesc(lh_2)
        hold on
    %     colorbar;
        plot(col,row,'r*');
        text(col, row, num2str(xlow + col*resolution) + "," + num2str(xlow + row*resolution));
        for i = [1,2,3,4]
            plot((P2(2*(i-1)+1)-xlow)/resolution,(P2(2*(i-1)+2)-xlow)/resolution,'bx');
        end
        view([90,-90]);
        ax = gca;
        ax.FontSize = 16; 
        xlabel("x1")
        ylabel("x2")
        xticks([1 20 40 60 80 ])
        xticklabels({xlow + 0, xlow+20*resolution, xlow+40*resolution, xlow+60*resolution, xlow+80*resolution})
        yticks([1 20 40 60 80 ])
        yticklabels({xlow + 0, xlow+20*resolution, xlow+40*resolution, xlow+60*resolution, xlow+80*resolution})
%         title("Configuration 2 and square loss with target at [0,0]")
        saveas(gcf, "figs/config_2_loss_analysis","epsc");
        saveas(gcf, "figs/config_2_loss_analysis","png");
        hold off
    end
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

    estimates = zeros(3,131);
    for time = [1:131]
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

    estimates = zeros(2,131);
    for time = [1:131]
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

% Load and plot
if 0
    load estimates_a.mat
    plot(estimates(1,:),estimates(2,:))
%     title("Position estimated with least squares grid search and sensor model 1")
    view([-90,90]);
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/position_grid_search_model_1","epsc");
    saveas(gcf, "figs/position_grid_search_model_1","png");
end
if 0
    load estimates_d.mat
    plot(estimates(1,:),estimates(2,:))
%     title("Position estimated with least squares grid search and sensor model 2")
    view([-90,90]);
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/position_grid_search_model_2","epsc");
    saveas(gcf, "figs/position_grid_search_model_2","png");
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
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/tracking_ekf_motion_1","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_1","png");
end

% Motion model 2
if 0
    z = sig(estimates');
    z.fs=2;
    S_location_1 = sensormod(@(t,x,u,th)[x(1,:); x(2,:)], [2,0,2,0]);
    S_location_1.pe = trace(reference_PE1.P)/3*eye(2);
    S_location_1.fs = 2;
    S_location_1_ctcv2d = addsensor(ctcv2d, S_location_1);
    x = ekf(S_location_1_ctcv2d, z);
    xplot2(x, 'conf', 90);
    view([-90 90])
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/tracking_ekf_motion_2","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_2","png");
end

% b) Reference TDOA sensor model
y_reference_tdoa = zeros(3,132);
for time =[1:132]
    y_reference_tdoa(:,time) = reference_tdoa(0,[estimates(1,time),estimates(2,time)],0,P1);
end

% Motion model 1
if 0
    z = sig(y_reference_tdoa', 2);
    z.t = z.t + 1;
    z.fs=2;
    z.Py = reference_PE1.P;
    S_reference_1_cv2d = addsensor(cv2d, S_reference_1);
    S_reference_1_cv2d.x0 = [1,1,1,1]';
    [x_ekf_reference_1_cv2d,V] = ekf(S_reference_1_cv2d, z);
    xplot2(x_ekf_reference_1_cv2d, 'conf', 90);
    view([-90 90])
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/tracking_ekf_motion_1_2","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_1_2","png");
end

% Motion model 2
if 0
    z = sig(y_reference_tdoa', 2);
    z.t = z.t + 1;
    z.fs=2;
    z.Py = reference_PE1.P;
    S_reference_1_ctcv2d = addsensor(ctcv2d, S_reference_1);
    S_reference_1_ctcv2d.x0 = [1,1,1,1,0]';
    [x_ekf_reference_1_ctcv2d,V] = ekf(S_reference_1_ctcv2d, z);
    xplot2(x_ekf_reference_1_ctcv2d, 'conf', 90);
    view([-90 90])
    ax = gca;
    ax.FontSize = 16; 
    xlabel("x1")
    ylabel("x2")
    saveas(gcf, "figs/tracking_ekf_motion_2_2","epsc");
    saveas(gcf, "figs/tracking_ekf_motion_2_2","png");
end
%% 7 Uncertainty
% How much can I perturb before really bad (maximum dist)
if 0
    for i = [1:20]
        i
        z = sig(y_reference_tdoa', 2);
        z.t = z.t + 1;
        z.fs=2;
        z.Py = reference_PE1.P;
        S_reference_1_perturbed = S_reference_1;
        S_reference_1_perturbed.th(1) = S_reference_1_perturbed.th(1) + i
        S_reference_1_ctcv2d = addsensor(ctcv2d, S_reference_1_perturbed);
        S_reference_1_ctcv2d.x0 = [1,1,1,1,0]';
        [x_ekf_reference_1_ctcv2d,V] = ekf(S_reference_1_ctcv2d, z);
        xplot2(x_ekf_reference_1_ctcv2d, 'conf', 90);
        view([-90 90])
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'xlabel',[])
        set(gca,'ylabel',[])
        saveas(gcf, ['figs/disturb_' num2str(i)],"epsc");
        saveas(gcf, ['figs/disturb_' num2str(i)],"png");
    end
end
