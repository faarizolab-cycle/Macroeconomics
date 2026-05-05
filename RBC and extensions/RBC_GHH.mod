% RBC Model with GHH Preferences
% Author: Francisco Arizola Blua
var C I K N Y W R A X;
varexo eps;

parameters beta delta alpha theta chi rho sigma sigma_eps;

% Calibration
beta  = 0.96;
delta = 0.025;
alpha = 0.36;
theta = 2.5;       % Calibrated to yield N_ss approx 0.33
chi   = 1;         % Inverse Frisch elasticity
rho   = 0.95;
sigma = 1;         % Logarithmic over the composite
sigma_eps = 0.01;

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
R_ss = 1/beta - 1 + delta;

% Exact SS for Labor (N) with GHH
N_ss = ( ((1 - alpha)/theta) * (alpha/R_ss)^(alpha/(1-alpha)) )^(1/chi);

% Ratios and Levels
K_N_ratio = (alpha / R_ss)^(1 / (1 - alpha));
K_ss = K_N_ratio * N_ss;
Y_ss = K_N_ratio^alpha * N_ss;
I_ss = delta * K_ss;
C_ss = Y_ss - I_ss;
W_ss = (1 - alpha) * Y_ss / N_ss;

% GHH Composite
X_ss = C_ss - (theta * N_ss^(1+chi)) / (1+chi);

model;
% 1. Composite Definition
X = C - (theta * N^(1+chi)) / (1+chi);

% 2. Euler equation (GHH)
X^(-sigma) = beta * X(+1)^(-sigma) * (R(+1) + 1 - delta);

% 3. Labor supply (GHH - No wealth effect)
theta * N^chi = W;

% 4. Capital accumulation
K = (1 - delta) * K(-1) + I;

% 5. Goods market
Y = C + I;

% 6. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% 7. Labor demand
W = (1 - alpha) * Y / N;

% 8. Capital demand
R = alpha * Y / K(-1);

% 9. TFP shock process
log(A) = rho * log(A(-1)) + eps;
end;

initval;
A = A_ss;
R = R_ss;
K = K_ss;
I = I_ss;
Y = Y_ss;
C = C_ss;
N = N_ss;
W = W_ss;
X = X_ss;
end;

steady;
check;

shocks;
var eps; stderr sigma_eps;
end;

stoch_simul(order=1, irf=40) A Y N C I W R;