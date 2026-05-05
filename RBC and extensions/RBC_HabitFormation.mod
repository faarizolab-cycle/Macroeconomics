% RBC Model with Internal Habit Formation
% Author: Francisco Arizola Blua 
var C I K N Y W R A Lambda;
varexo eps;

parameters beta delta alpha theta phi rho sigma;

% Calibration
beta  = 0.96;
delta = 0.025;
alpha = 0.36;
theta = 1;     % Preference for leisure
rho   = 0.95;
sigma = 0.01;

% --- MACRO-PROCESADOR DE DYNARE ---
% Si MATLAB no envía un valor, se asume 0.9 por defecto.
@#ifndef phi_val
    @#define phi_val = 0.9
@#endif
phi = @{phi_val};
% ----------------------------------

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
R_ss = 1/beta - 1 + delta;

% Ratios
K_N_ratio = (alpha * A_ss / R_ss)^(1 / (1 - alpha));
Y_N_ratio = A_ss * K_N_ratio^alpha;
W_ss      = (1 - alpha) * Y_N_ratio;
I_N_ratio = delta * K_N_ratio;
C_N_ratio = Y_N_ratio - I_N_ratio;

% Exact SS for N with Internal Habit:
term_X = (W_ss * (1 - beta*phi)) / (theta * C_N_ratio * (1 - phi));
N_ss   = term_X / (1 + term_X);

% Niveles
K_ss = K_N_ratio * N_ss;
Y_ss = Y_N_ratio * N_ss;
I_ss = I_N_ratio * N_ss;
C_ss = C_N_ratio * N_ss;

% Lambda SS
Lambda_ss = (1 - beta*phi) / (C_ss * (1 - phi));

model;
% 1. Marginal Utility of Consumption (Lambda)
Lambda = (1 / (C - phi*C(-1))) - beta*phi*(1 / (C(+1) - phi*C));

% 2. Euler equation
Lambda = beta * Lambda(+1) * (R(+1) + 1 - delta);

% 3. Labor supply 
theta / (1 - N) = Lambda * W;

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
Lambda = Lambda_ss;
end;

steady;
check;

shocks;
var eps; stderr sigma;
end;

stoch_simul(order=1, irf=40) A Y N C I W R;