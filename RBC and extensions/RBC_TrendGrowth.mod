% RBC Model with Stochastic Trend Growth
% Author: Francisco Arizola Blua
var C_hat I_hat K_hat N Y_hat W_hat R lambda z g;
varexo eps_z eps_g;

parameters beta delta alpha psi sigma rho rho_g mu_g sigma_z sigma_g;

% Calibration
beta   = 0.96;
delta  = 0.025;
alpha  = 0.64;     % Labor share (since Y = K^(1-alpha) * N^alpha)
psi    = 0.33;     % Weight on consumption in utility
sigma  = 2;        % Coefficient of relative risk aversion
rho    = 0.95;     % Persistence of stationary TFP shock
mu_g   = 1.005;    % Steady state gross growth rate (e.g., 0.5% per quarter)
sigma_z = 0.01;
sigma_g = 0.005;

% --- MACRO-PROCESADOR DE DYNARE ---
% Recibe el valor de rho_g desde MATLAB (persistencia del shock de crecimiento).
@#ifndef rhog_val
    @#define rhog_val = 0.5
@#endif
rho_g = @{rhog_val};
% ----------------------------------

% --------------------------------------------------------
% Steady State computation (Detrended variables)
% --------------------------------------------------------
z_ss = 0;
g_ss = mu_g;

% Euler equation SS solved for R
R_ss = (mu_g^(1 - psi*(1-sigma))) / beta - 1 + delta;

% Ratios
K_N_ratio = ((1 - alpha) / R_ss)^(1 / alpha);
Y_N_ratio = K_N_ratio^(1 - alpha);
W_hat_ss  = alpha * Y_N_ratio;
I_N_ratio = (mu_g - (1 - delta)) * K_N_ratio;
C_N_ratio = Y_N_ratio - I_N_ratio;

% Exact SS for Labor (N)
term_denom = 1 - (mu_g - 1 + delta) * ((1 - alpha) / R_ss);
term_num   = (psi * alpha) / (1 - psi);
N_ss       = term_num / (term_denom + term_num);

% Levels
K_hat_ss = K_N_ratio * N_ss;
Y_hat_ss = Y_N_ratio * N_ss;
I_hat_ss = I_N_ratio * N_ss;
C_hat_ss = C_N_ratio * N_ss;

% Lambda SS
lambda_ss = psi * C_hat_ss^(psi*(1-sigma)-1) * (1 - N_ss)^((1-psi)*(1-sigma));

model;
% 1. Marginal Utility of Consumption (lambda)
lambda = psi * C_hat^(psi*(1-sigma)-1) * (1 - N)^((1-psi)*(1-sigma));

% 2. Euler equation (Detrended)
% Note the exact scaling factor for future marginal utility
lambda = beta * lambda(+1) * g(+1)^(psi*(1-sigma)-1) * (R(+1) + 1 - delta);

% 3. Intratemporal Labor Supply
((1-psi)/psi) * (C_hat / (1-N)) = W_hat;

% 4. Production function (Detrended)
Y_hat = exp(z) * K_hat(-1)^(1-alpha) * N^alpha;

% 5. Labor demand
W_hat = alpha * Y_hat / N;

% 6. Capital demand
R = (1-alpha) * Y_hat / K_hat(-1);

% 7. Capital accumulation (Detrended law of motion)
K_hat = (1/g) * ((1-delta)*K_hat(-1) + I_hat);

% 8. Goods market
Y_hat = C_hat + I_hat;

% 9. Stationary TFP shock
z = rho * z(-1) + eps_z;

% 10. Trend Growth shock
log(g) = (1-rho_g)*log(mu_g) + rho_g*log(g(-1)) + eps_g;
end;

initval;
z = z_ss;
g = g_ss;
R = R_ss;
K_hat = K_hat_ss;
I_hat = I_hat_ss;
Y_hat = Y_hat_ss;
C_hat = C_hat_ss;
N = N_ss;
W_hat = W_hat_ss;
lambda = lambda_ss;
end;

steady;
check;

shocks;
var eps_z; stderr sigma_z;
var eps_g; stderr sigma_g;
end;

% Generamos las IRFs
stoch_simul(order=1, irf=40) Y_hat C_hat I_hat N W_hat R K_hat g;