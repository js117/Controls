% KalmanTestRoutine
% This function simulates a Kalman filter smoothing out y_data, and
% compares results graphically of original y_data vs. y_data_smoothed.


function [x_kalman, y_error] = KalmanTestRoutine(S0, x0bar, W, wbar, V, vbar, A, C, y_data, numTimeSteps)

    % need some dimensions.. function users should ensure their problem
    % data corresponds to the below conventions.
    n = size(A,1); %A is nxn
    p = size(C,1); %C is pxn
    
    x_kalman = zeros(n, numTimeSteps); %i.e. every column will be an updated state snapshot
    y_error = zeros(p, numTimeSteps); %ditto for y; note this is also the size of y_data as given
    
    x_hat_previous = x0bar;
    x_hat_measured = zeros(size(x0bar));
    x_hat_updated = zeros(size(x0bar));
    S_previous = S0;
    S_measured = zeros(size(S0));
    S_updated = zeros(size(S0));
    
    % Simulation. For a real-time system, this loop runs continually.
    for i=1:numTimeSteps
  
       % Apply measurement and time updates: 
       [x_hat_measured, S_measured] = KalmanMeasurementUpdate(x_hat_previous, S_previous, C, V, y_data(:,i), vbar);
       [x_hat_updated, S_updated] = KalmanTimeUpdate(x_hat_measured, S_measured, A, W, wbar);
       
       % Recrod the Kalman filtered state estimate, and the measurement error: 
       % (in application code, this is where we look at quantities of interest)
       x_kalman(:,i) = x_hat_updated; %or take a look at x_hat_updated
       y_error(:,i) = y_data(:,i) - C*x_hat_measured - vbar; 
       
       % Update recursively for next iteration:
       S_previous = S_updated;
       x_hat_previous = x_hat_updated;
       
    end
    
end