% Riccati recursion script
% create A,B,C,Q,Qf,R and set N in workspace
% set as well x0 to compute the full state and input sequence

n = size(A,1); % A is size n x n
m = size(B,2); % B is size n x m

Pmats = zeros(n,n,N+1); % storing each iteration
Kmats = zeros(m,n,N); % storing each iteration
Pmats(:,:,N+1) = Qf; % final condition
convergenceTest = zeros(N,3); %assuming m at least 3
u_lqr = zeros(N,m);
x_opt = zeros(N+1,n);
if (exist('x0','var'))
   x_opt(1,:) = x0; 
end

for t = N+1:-1:2
   Pt = Pmats(:,:,t);
   term3 = (inv(R + B'*Pt*B))*B'*Pt*A; %this term is repeated in P(t-1), Kt
   Pmats(:,:,t-1) = Q + A'*Pt*A - A'*Pt*B*term3;
   Kmats(:,:,t-1) = -term3;
   
   [U,S,V] = svd(-term3); % check convergence via singular values
   
   %convergenceTest(t-1, 1) = S(1,1); %what to check in general
   %below: for ee363 pset 1 #1
   convergenceTest(t-1,1) = -term3(1);
   convergenceTest(t-1,2) = -term3(2);
   convergenceTest(t-1,3) = -term3(3);
end

u_lqr(1,:) = Kmats(:,:,1)*x0;
for t = 1:N %loop forward to calculate states
   if (exist('x0','var'))
       u_lqr(t,:) = Kmats(:,:,t)*x_opt(t,:)'; %set next u(t) val
       x_opt(t+1,:) = A*x_opt(t,:)' + B*u_lqr(t,:)';
   end 
end

% plot the largest singular value of K to see how much it changes:

% for pset, plot 3 values of Kt to see when we converge
figure; plot(convergenceTest(:,1));
figure; plot(convergenceTest(:,2));
figure; plot(convergenceTest(:,3));

figure; plot(x_opt(:,1));
figure; plot(x_opt(:,2));
figure; plot(x_opt(:,3));
figure; plot(u_lqr(:,1));