% RBC Model with Indivisible Labor (Hansen-Rogerson)
% Author: Francisco Arizola Blua
var C K N Y W R A I;
varexo eps;

parameters beta alpha delta rho sigma B;

% Calibration
beta   = 0.99;     % Discount factor
alpha  = 0.33;     % Capital share
delta  = 0.025;    % Depreciation rate
rho    = 0.95;     % TFP persistence
sigma  = 0.01;     % Std. dev. of TFP shock
B      = 2.63;     % Disutility of labor

% Steady state computation (outside Dynare block)
A_ss = 1;
R_ss = 1/beta - 1 + delta;
K_N_ratio = (alpha * A_ss / R_ss)^(1 / (1 - alpha));
Y_N_ratio = A_ss * K_N_ratio^alpha;
W_ss = (1 - alpha) * Y_N_ratio;
C_N_ratio = (1 - alpha)/B * Y_N_ratio;

N_ss = Y_N_ratio / (C_N_ratio + delta * K_N_ratio);
K_ss = K_N_ratio * N_ss;
Y_ss = Y_N_ratio * N_ss;
C_ss = C_N_ratio * N_ss;
I_ss = Y_ss - C_ss;
W_ss = W_ss;
R_ss = R_ss;

% Model equations
model;
% 1. Euler equation
1/C = beta * (1/C(+1)) * (1 + R(+1) - delta);
% 2. Labor supply (indivisible)
B = (1/C) * W;
% 3. Capital accumulation
K = (1 - delta) * K(-1) + I;
% 4. Goods market identity (defines Investment)
I = Y - C;
% 5. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);
% 6. Labor demand
W = (1 - alpha) * Y / N;
% 7. Capital demand
R = alpha * Y / K(-1);
% 8. TFP process
log(A) = rho * log(A(-1)) + eps;
end;

initval;
A = A_ss;
R = R_ss;
K = K_ss;
Y = Y_ss;
N = N_ss;
C = C_ss;
I = I_ss;
W = W_ss;
end;

steady;
check;

shocks;
var eps; stderr sigma;
end;

stoch_simul(order=1,irf=40) A Y N C I W R;