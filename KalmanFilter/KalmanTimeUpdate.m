% -------------------
% --- Time Update ---
% ----------------------------
% x_hat(t+1|t) = A*x_hat(t|t)
% 
% S(t+1|t) = A*S(t|t)*A' + W
% ----------------------------
% Extra modification: including wbar as the contribution due to mean(w)!=0

function [x_hat_updated, S_updated] = KalmanTimeUpdate(x_hat_measured, S_measured, A, W, wbar)

    x_hat_updated = A*x_hat_measured + wbar;
    
    S_updated = A*S_measured*A' + W; %doesn't change due to wbar

end