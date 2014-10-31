% --------------------------
% --- Measurement Update ---
% --------------------------------------------------------------------------------------
% x_hat(t|t) = x_hat(t|t-1) + S(t|t-1)*C'*(C*S(t|t-1)*C'+V).inv*(y(t) - C*x_hat(t|t-1))
% 
% S(t|t) = S(t|t-1) - S(t|t-1)*C'*(C*S(t|t-1)*C'+V).inv*C*S(t|t-1)
% --------------------------------------------------------------------------------------
% Extra modification: including vbar as the contribution due to mean(v)!=0

function [x_hat_measured, S_measured] = KalmanMeasurementUpdate(x_hat_previous, S_previous, C, V, yt, vbar)

    tempInv = inv(C*S_previous*C'+V);

    x_hat_measured = x_hat_previous + S_previous*C'*tempInv*(yt - C*x_hat_previous - vbar); %vs just (y(t) - C*x_hat(t|t-1))

    S_measured = S_previous - S_previous*C'*tempInv*C*S_previous;
    
end