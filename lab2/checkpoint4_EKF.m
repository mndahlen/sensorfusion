function [xhat, meas] = checkpoint2_EKF(fname, calAcc, calGyr, calMag)

  %% Setup necessary infrastructure
  import('se.hendeby.sensordata.*');  % Used to receive data.

  DISPLAY_FREQ = 10; % [Hz]  Frequency of update of the visualization

  %% Filter settings
  t0 = [];  % Initial time (initialize on first data received)
  nx = 4;
  % Add your filter settings here.
  % Current filter state.
  x = [1; 0; 0; 0];
  P = eye(nx, nx);
  prev_gyr = [0,0,0]';
  prev_acc = [0,0,-calAcc.m(3)]';
  prev_mag = calMag.m;

  % Saved filter states.
  xhat = struct('t', zeros(1, 0),...
                'x', zeros(nx, 0),...
                'P', zeros(nx, nx, 0));

  meas = struct('t', zeros(1, 0),...
                'acc', zeros(3, 0),...
                'gyr', zeros(3, 0),...
                'mag', zeros(3, 0),...
                'orient', zeros(4, 0));
  try
    %% Create data link
    err_hint = ['Unsuccessful loading of sensordata.jar!\n',...
        'Make sure to run startup.m before attempting to run this function.'];
    if nargin == 0 || isempty(fname)
      server = StreamSensorDataReader(3400);
      err_hint = ['Unsuccessful connecting to client!\n',...
        'Make sure to start streaming from the phone *after* starting this function.'];
    else
      server = FileSensorDataReader(fname);
      err_hint = sprintf(['Unsuccessful reading data from file!\n',...
        'Make sure ''%s'' is the correct path to the log file.'], fname);
    end
    % Makes sure to resources are returned.
    sentinel = onCleanup(@() server.stop());

    server.start();  % Start data reception.
    clear err_hint
  catch e
    e = e.addCause(MException(sprintf('%s:UnableToConnect', mfilename), err_hint));
    rethrow(e);
  end

  % Used for visualization.
  figure(1);
  subplot(1, 2, 1);
  ownView = OrientationView('Own filter', gca);  % Used for visualization.
  ownView.activateKeyboardCallback;
  googleView = [];
  tnextdisp = 0;  % Next time to update the visualization

  %% Filter loop
  % Repeat while data is available and q hasn't been pressed
  while server.status() && ~ownView.quit
    % Get the next measurement set, assume all measurements
    % within the next 5 ms are concurrent (suitable for sampling
    % in 100Hz).
    data = server.getNext(5);

    if isnan(data(1))  
      % No new data received
      continue;
    end
    % Extract current time
    t = data(1)/1000;  

    % Initialize t0
    if isempty(t0)  
      t0 = t;
    end
    
    % Gyroscope 
    gyr = data(1, 5:7)';
    if ~any(isnan(gyr)) 
        % Gyro measurements are available.
        gyr;
    else
        gyr = prev_gyr;
    end
    prev_gyr = gyr;
    
    % Accelerometer
    acc = data(1, 2:4)';
    if ~any(isnan(acc))  
        % Acc measurements are available.
        % Outlier Rejection (Assuming dominated by g)
        gaussian_exponent = (acc(3)-calAcc.m(3))^2/calAcc.R(3,3);
        if gaussian_exponent > 5e5
            ownView.setAccDist(1);
        else
            ownView.setAccDist(0)
            [x, P] = mu_g(x, P, acc, calAcc.R, [0,0,calAcc.m(3)]');
        end
        % Measurement Update
    else
        acc = prev_acc;
    end
    prev_acc = acc;
    
    % Magnetometer
    mag = data(1, 8:10)';
    if ~any(isnan(mag))
        % Mag measurements are available.
        % Outlier Rejection (Assuming dominated by g)
        gaussian_exponent = (mag-calMag.m)'*inv(calMag.R)*(mag-calMag.m);
        if gaussian_exponent > 1e7
            ownView.setMagDist(1);
        else
            [x, P] = mu_m(x, P, mag, calMag.R, [0,sqrt(calMag.m(1)^2+calMag.m(2)^2),calMag.m(3)]');
            ownView.setMagDist(0)
        end
        % Measurement Update
    else
        mag = prev_mag;
    end
    prev_mag = mag;
    
    % Time update
    [x, P] = tu_qw_1(x, P, gyr, 0.01, calGyr.R);
    % Normalize
    [x, P] = mu_normalizeQ(x, P);

    % Google's orientation estimate.
    orientation = data(1, 18:21)';  
   
    % Visualize result
    if t >= tnextdisp
      tnextdisp = t + 1/DISPLAY_FREQ;  % Next vizualization update
      setOrientation(ownView, x(1:4));
      title(ownView, 'OWN', 'FontSize', 16);
      if ~any(isnan(orientation))
        if isempty(googleView)
          subplot(1, 2, 2);
          % Used for visualization.
          googleView = OrientationView('Google filter', gca);
        end
        setOrientation(googleView, orientation);
        title(googleView, 'GOOGLE', 'FontSize', 16);
      end
    end

    % Save estimates
    xhat.x(:, end+1) = x;
    xhat.P(:, :, end+1) = P;
    xhat.t(end+1) = t - t0;
    meas.t(end+1) = t - t0;
    % Save measurements
    meas.acc(:, end+1) = acc;
    meas.gyr(:, end+1) = gyr;
    meas.mag(:, end+1) = mag;
    meas.orient(:, end+1) = orientation;
  end
end