% % KalmanTestScript
% 
% In this toy example, we will consider applying a Kalman filter to a sensor
% problem: the determination of position from accelerometer readings, in one axis.
%
% Let the angle measurement we wish to estimate be x(t); let the gyroscope
% reading be g(t); let the time between samples be dt. Then we have:
%
% [x(t+1), gp(t+1)]' = [1, dt; 0, 1]*[x(t) g(t)]' + w(t)
% (we have a prediction of next angular velocity gp that we don't really use)
% y(t) = [0, 1]*[x(t) g(t)]' + v(t)
%
% Where w(t) and v(t) are process and measurement noise, respectively. We
% will simulate a very small process noise (representing noise between
% system time steps) and the measurement noise being the noise of the gyro.
%
% 
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. create a true g(t) vector:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numTimeSteps = 1000; %e.g. at sample rate 100Hz, this is 1 second worth of data
randValues = randn(numTimeSteps,1); %mean 0, std 1
g = zeros(numTimeSteps, 1);
gtrue = zeros(numTimeSteps, 1);
gsigma = 150;
gUniNoise = 5;
gbar = 100 + gUniNoise/2; %add some mean value offset
% We add gUniNoise/2 b/c rand() returns a val from (0,1) so on average this
% produces a mean offset of half the range.

% Create sample data split into thirds: rise, steady, fall back down.
% For example this could be an actuator performing a repetitive task, and
% the gyro measures the motor's applied angular velocity.
% Do this twice. 

% Note: randValues are Gaussian which fits the Kalman measurement model
% well, so even at large gsigma we can filter quite well. OTOH, the noise
% from MATLAB's rand() is uniform, and this particular noise is harder to
% deal with and causes more error; this makes sense as the KF is derived
% with respect to Gaussian noise. One could imagine re-deriving it for
% different noise models (would be a good TODO).
for i=1:floor(numTimeSteps/6)
   gtrue(i+1) = gtrue(i) + 5; %rising linearly
   g(i+1) = gtrue(i+1) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1);  
end
for i=ceil(numTimeSteps/6):floor(2*numTimeSteps/6)
   gtrue(i+1) = 0; %stall, i.e. angular velocity of 0
   g(i+1) = gtrue(i) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1);
end
for i=ceil(2*numTimeSteps/6):floor(3*numTimeSteps/6)
   gtrue(i+1) = gtrue(i) - 5; %negative angular velocity over last third to get back to start
   g(i+1) = gtrue(i+1) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1); 
end
gtrue(ceil(3*numTimeSteps/6)) = 0; % remember angular velocity is directional, doesn't start where left off
for i=ceil(3*numTimeSteps/6):floor(4*numTimeSteps/6)
   gtrue(i+1) = gtrue(i) + 5; %rising linearly
   g(i+1) = gtrue(i+1) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1);  
end
for i=ceil(4*numTimeSteps/6):floor(5*numTimeSteps/6)
   gtrue(i+1) = 0; %stall, i.e. angular velocity of 0
   g(i+1) = gtrue(i) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1);
end
for i=ceil(5*numTimeSteps/6):numTimeSteps-1
   gtrue(i+1) = gtrue(i) - 5; %negative angular velocity over last third to get back to start
   g(i+1) = gtrue(i+1) + gbar + gsigma*randValues(i) + gUniNoise*rand(1,1); 
end
xaxis = (1:numTimeSteps)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Simulate the true state model (i.e. without noise:)
% [x(t+1), gp(t+1)]' = [1, dt; 0, 1]*[x(t) g(t)]' 
%
% y(t) = [0, 1]*[x(t) g(t)]' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 1e-2; %e.g. 10ms between samples; freq = 100Hz
A = [1, dt; 0, 1];
C = [0, 1];
x0 = [0, 0]';
n = 2;
p = 1;
x_state_true = zeros(n, numTimeSteps); %every column is a state x(t)
x_state_noisy = zeros(n, numTimeSteps);
y_measured_true = zeros(p, numTimeSteps); %every column is a state x(t)
y_measured_noisy = zeros(p, numTimeSteps);

for t=1:numTimeSteps-1
   % true reference
   x_state_true(:,t+1) = A*[x_state_true(1,t), gtrue(t)]';
   y_measured_true(:,t) = C*[x_state_true(1,t), gtrue(t)]';
   % with noise model
   x_state_noisy(:,t+1) = A*[x_state_noisy(1,t), g(t)]';
   y_measured_noisy(:,t) = C*[x_state_noisy(1,t), g(t)]';
end
y_measured_true(:,numTimeSteps) = C*[x_state_true(1,numTimeSteps), gtrue(numTimeSteps)]';

angles_true = x_state_true(1,:);
angles_noisy = x_state_noisy(1,:);

% Observe that nasty integration error!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Kalman Filter Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0bar = x0;
S0 = 0; % we start (and expect to start) from a state of rest for this process
V = 0.001*gsigma; % in practice we don't know gsigma exactly, so let's say we underestimate it
vbar = gbar;
W = 0.001*[1,0;0,10]; %"process noise"
wbar = 0;
y_data = zeros(p,numTimeSteps);
y_data(1,:) = g;
%[x_kalman, y_error] = KalmanTestRoutine(S0, x0bar, W, wbar, V, vbar, A, C, y_data, numTimeSteps);

% Below code adapted from KalmanTestRoutine:
x_kalman = zeros(n, numTimeSteps); %i.e. every column will be an updated state snapshot
y_error = zeros(p, numTimeSteps); %ditto for y; note this is also the size of y_data as given
x_hat_previous = x0bar;
x_hat_measured = zeros(size(x0bar));
x_hat_updated = zeros(size(x0bar));
S_previous = S0;
% Simulation. For a real-time system, this loop runs continually.
for i=1:numTimeSteps

   % Apply measurement and time updates: 
   [x_hat_measured, S_measured] = KalmanMeasurementUpdate(x_hat_previous, S_previous, C, V, y_data(1,i), vbar);
   [x_hat_updated, S_updated] = KalmanTimeUpdate(x_hat_measured, S_measured, A, W, wbar);

   % Recrod the Kalman filtered state estimate, and the measurement error: 
   % (in application code, this is where we look at quantities of interest)
   x_kalman(:,i) = x_hat_updated; %or take a look at x_hat_updated
   y_error(:,i) = y_data(:,i) - C*x_hat_measured - vbar; 

   % Update recursively for next iteration:
   S_previous = S_updated;
   x_hat_previous = x_hat_updated;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Plot and review results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the velocities: true, simulated noisy, Kalman filtered:
figure; plot(xaxis, gtrue, 'k', xaxis, g, 'r', xaxis, x_kalman(2,:), 'g');

% Plot and look at the actual + noisy angle readings:
figure; plot(xaxis, angles_true, 'k', xaxis, angles_noisy, 'r', xaxis, x_kalman(1,:), 'g');

