%% TSRT14, Lab 2: Orientation
% _This lab report describes the experiences made while working with lab 2
% in TSRT14 Sensor Fusion._
%
%% *Group members:*
%
% # Martin Dahl (981223-2677), marda545
%
%% _Usage_
% _This file is intended to be a template for reporting the results on lab

% if inpublish  % Load saved data when publishing.
%   load DATAFILE
% end

%% 2. Get to know your data
% The accelerometer has a component from the earth gravitational field
% which is constantly around 9.82.
% The magnetometer reacts somewhat to movement and almost always has non
% zero values
% The gyroscope measurements are 0 if the phone is not moving. 

% if ~inpublish  % Don't recollect data during publish
%   [xhat2, meas2] = ekfFilter();
%   save DATAFILE -append xhat2 meas2
% end
%%
%%
% *Result*
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


%% 3. Add the EKF time update step
% _*Include your time update function* from the
% preparations here:
function [x, P] = tu_qw(x, P, omega, T, Rw)
    S_q = Sq(x);
    S_omega = Somega(omega);
    I = eye(size(S_omega));

    F = (I + 0.5*S_omega*T);
    G = 0.5*T*S_q;
    Q = Rw;
    
    x = F*x;
    P = F*P*(F') + G*Q*(G');
end
%
% <include>tu_qw.m</include>
%
%%
% Time update is based on measurement, so we use the measurement noise
% covariance.
%
%%
% *Result*
%
% _Run the indicated code below to generate results to plot._
% if ~inpublish  % Don't recollect data during publish
%   [xhat3, meas3] = ekfFilter();
%   save DATAFILE -append xhat3 meas3
% end
%%
if 0
    figure
    visDiff(xhat3, meas3);
end

%%
% *_Shortly describe your observations:_*
% The phone reacts very well to time update using gyroscope measurements.
% We now measure angular velocity which can be used for time update.
% Furthermore, the gyroscope helps the phone know if it is static or not.
% When we start flat we start from the calibration point. When we don't we
% start from something unexpected. We have not linearized around this
% point.

%% 4. Add the EKF accelerometer measurement update step
function [x, P] = mu_g(x, P, yacc, Ra, g0)
    Q = Qq(x);
    [dQ0, dQ1, dQ2, dQ3] = dQqdq(x);
    H = (Q')*g0; 
    dH = [dQ0'*g0, dQ1'*g0, dQ2'*g0, dQ3'*g0]; % Seems to drift?

    S = Ra + dH*P*(dH');
    K = P*(dH')*S^(-1);
    eps = yacc - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end
%
% <include>mu_g.m</include>
%
%%
% We include the acceleration measurement covariance and g0 (mean of
% acceleration where we assumate gravitation dominates).

%% Note: USe [0,0,g0] and subtract the estimated error in component 1 and 2 form measurements
%
%%
% *Result*
% if ~inpublish  % Don't recollect data during publish
%   [xhat4, meas4] = ekfFilter();
%   save DATAFILE -append xhat4 meas4
% end
%%
if 0
    figure
    visDiff(xhat4, meas4);
end
%%
% Ability to update measurement with acceleration, for which we can use
% gravity and this is useful if the phone gets lost in orientation. If we
% shake it heavily we deviate from the point of linearization too quickly. 
% It's difficult to observe the effect of sliding the phone on a surface.

%% 5. Add accelerometer outlier rejection
% We assume accelerations follow a gaussian. If the gaussian PDF exponent
% value is above a threshold we say its an outlier. If we find an outlier
% we simply ignore the measurment and don't perform the measurement update.
%%
if 0
    gaussian_exponent = (acc(3)-calAcc.m(3))^2/calAcc.R(3,3);
    if gaussian_exponent > 5e5
        ownView.setAccDist(1);
    else
        ownView.setAccDist(0)
        [x, P] = mu_g(x, P, acc, calAcc.R, [0,0,calAcc.m(3)]');
    end
end
%%
% *Result*
%
% _Run the indicated code below to generate results to plot._
if ~inpublish  % Don't recollect data during publish
  [xhat5, meas5] = ekfFilter();
  save DATAFILE -append xhat5 meas5
end
%%
if 0
    figure
    visDiff(xhat5, meas5);
end
%%
% _*Shortly describe your observations:*_
%
% * _What happens when you shake or quickly slide the phone on surface?
%   Why?_
%% NOTE: Should be an outlier since we are not expected external forces
% 

%% 6. Add the EKF magnetometer measurement update step
% _*Include your magnetometer measurement update function* from the
% preparations here:_
function [x, P] = mu_m(x, P, ymag, Rm, m0)
    Q = Qq(x);
    [dQ0, dQ1, dQ2, dQ3] = dQqdq(x);
    H = (Q')*m0; 
    dH = [dQ0'*m0, dQ1'*m0, dQ2'*m0, dQ3'*m0];

    S = Rm + dH*P*(dH');
    K = P*(dH')*S^(-1); % (S\eye(length(S))
    eps = ymag - H;

    x = x + K*eps;
    P = P - P*(dH')*S^(-1)*dH*P;
end

% <include>mu_m.m</include>
%
%%
% _*Motivate parameter choices:*_
% Like before
%
%
%%
% *Result*
%
% _Run the indicated code below to generate results to plot._
if ~inpublish  % Don't recollect data during publish
  [xhat6, meas6] = ekfFilter();
  save DATAFILE -append xhat6 meas6
end
%%
if 0 
    figure
    visDiff(xhat6, meas6);
end
%%
% _*Shortly describe your observations:*_
% We can now get back to our actual orientation since the earth magnetic
% field is static, however, its not as good as using gravitation. 
% But the magnetic field is in another axis than gravitation. If the
% phone is clsoe to a magnet it flips out, so we need some outlier
% rejection.
%% Note: Compared to accelerometer with g wwe now have constant component in all axes

%% 7. Add magnetometer outlier rejection
% Similar to accelerometer outlier rejection we assume a gaussian
% distribution and check for too high PDF exponent values. We use a
% threshold to determine if something is an outlier. If we find something
% to be an outlier we ignore the measurement and don't perform the
% measurement update using the measurement.
%%
if 0
    gaussian_exponent = (mag-calMag.m)'*inv(calMag.R)*(mag-calMag.m);
    if gaussian_exponent > 1e7
        ownView.setMagDist(1);
    else
        [x, P] = mu_m(x, P, mag, calMag.R, [0,sqrt(calMag.m(1)^2+calMag.m(2)^2),calMag.m(3)]');
        ownView.setMagDist(0)
    end
end
%
%%
% *Result*
%
% if ~inpublish  % Don't recollect data during publish
%   [xhat7, meas7] = ekfFilter();
%   save DATAFILE -append xhat7 meas7
% end
%%
if 0
    figure
    visDiff(xhat7, meas7);
end
%%
% When we have outlier rejection it no longer flips out close a magnet. It
% is more stable.

%% 8. Test your filter without gyroscope measurements
%%
% *Result*
%
% if ~inpublish  % Don't recollect data during publish
%   [xhat8, meas8] = ekfFilter();
%   save DATAFILE -append xhat8 meas8
% end
if 0
    figure
    visDiff(xhat8, meas8);
end
%%
% It moves very slowly and does not react strongly to the phone rotating.
% The orientation now has to be found using only measurement updates which
% takes more time. If the model for other measurements are better then we
% could potentially use the time update on another measurement.
%% Note: Increase process noise to increase sensitivity using only MU.


%% APPENDIX: Main loop
%
% <include>ekfFilter.m</include>
%
