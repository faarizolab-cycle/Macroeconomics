% RBC Model with KPR Preferences
% Author: Francisco Arizola Blua
var C I K N Y W R A;
varexo eps;

parameters beta delta alpha psi rho sigma sigma_eps;

% Calibration
beta  = 0.96;
delta = 0.025;
alpha = 0.36;
psi   = 2;
rho   = 0.95;
sigma_eps = 0.01;

% --- MACRO-PROCESADOR DE DYNARE ---
% Recibe el valor de sigma desde MATLAB. Por defecto es 1 (Log).
@#ifndef sigma_val
    @#define sigma_val = 1
@#endif
sigma = @{sigma_val};
% ----------------------------------

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
R_ss = 1/beta - 1 + delta;

% Exact SS for Labor (N) with KPR
N_ss = (1 - alpha) / (psi * (1 - delta * alpha / R_ss) + (1 - alpha));

% Ratios and Levels
K_N_ratio = (alpha / R_ss)^(1 / (1 - alpha));
K_ss = K_N_ratio * N_ss;
Y_ss = K_N_ratio^alpha * N_ss;
I_ss = delta * K_ss;
C_ss = Y_ss - I_ss;
W_ss = (1 - alpha) * Y_ss / N_ss;

model;
% 1. Euler equation (KPR)
(C*(1-N)^psi)^(-sigma) * (1-N)^psi = beta * (C(+1)*(1-N(+1))^psi)^(-sigma) * (1-N(+1))^psi * (R(+1) + 1 - delta);

% 2. Labor supply (KPR)
(psi * C) / (1 - N) = W;

% 3. Capital accumulation
K = (1 - delta) * K(-1) + I;

% 4. Goods market
Y = C + I;

% 5. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% 6. Labor demand
W = (1 - alpha) * Y / N;

% 7. Capital demand
R = alpha * Y / K(-1);

% 8. TFP shock process
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
end;

steady;
check;

shocks;
var eps; stderr sigma_eps;
end;

stoch_simul(order=1, irf=40) A Y N C I W R;