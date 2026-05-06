% RBC Model with Money in the Utility Function (MIU)
% Author: Francisco Arizola Blua
var C I K N Y W R m i pi mu A;
varexo eps_a eps_mu;

parameters beta delta alpha theta chi phi rho_a rho_mu mu_bar sigma_a sigma_mu;

% Calibration
beta     = 0.96;
delta    = 0.025;
alpha    = 0.36;
theta    = 1;
chi      = 1;
phi      = 0.1;       % Weight on real money balances in utility
rho_a    = 0.95;
rho_mu   = 0.50;      % Persistence of money growth shock
mu_bar   = 1.01;      % Steady-state gross money growth rate (~1% inflation)
sigma_a  = 0.01;
sigma_mu = 0.005;

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss  = 1;
mu_ss = mu_bar;
pi_ss = mu_bar - 1;

R_ss = 1/beta - 1 + delta;
i_ss = mu_bar/beta - 1;

% Ratios
K_N_ratio = (alpha / R_ss)^(1 / (1 - alpha));
Y_N_ratio = K_N_ratio^alpha;
W_ss      = (1 - alpha) * Y_N_ratio;
I_N_ratio = delta * K_N_ratio;
C_N_ratio = Y_N_ratio - I_N_ratio;

% Exact SS for Labor (N)
N_ss = ( (1 - alpha) / (theta * (1 - delta * alpha / R_ss)) )^(1 / (1 + chi));

% Levels
K_ss = K_N_ratio * N_ss;
Y_ss = Y_N_ratio * N_ss;
I_ss = I_N_ratio * N_ss;
C_ss = C_N_ratio * N_ss;

% Real money balances
m_ss = phi * ((1 + i_ss) * C_ss) / i_ss;

model;
% 1. Euler equation (Capital)
1/C = beta * (1/C(+1)) * (R(+1) + 1 - delta);

% 2. Labor supply
theta * N^chi = W / C;

% 3. Money demand
m = phi * ((1 + i) * C) / i;

% 4. Nominal Euler Equation (Fisher relation)
1/(1 + i) = beta * (C / (C(+1) * (1 + pi(+1))));

% 5. Capital accumulation
K = (1 - delta) * K(-1) + I;

% 6. Resource constraint
Y = C + I;

% 7. Production function and factor prices
Y = A * K(-1)^alpha * N^(1 - alpha);

% 8. Labor demand
W = (1 - alpha) * Y / N;

% 9. Capital demand
R = alpha * Y / K(-1);

% 10. Money market clearing (in real terms)
m * (1 + pi) = mu * m(-1);

% 11. TFP shock process
log(A) = rho_a * log(A(-1)) + eps_a;

% 12. Money growth process
log(mu) = (1 - rho_mu)*log(mu_bar) + rho_mu*log(mu(-1)) + eps_mu;
end;

initval;
A  = A_ss;
mu = mu_ss;
pi = pi_ss;
i  = i_ss;
R  = R_ss;
K  = K_ss;
I  = I_ss;
Y  = Y_ss;
C  = C_ss;
N  = N_ss;
W  = W_ss;
m  = m_ss;
end;

steady;
check;

shocks;
var eps_a;  stderr sigma_a;
var eps_mu; stderr sigma_mu;
end;

stoch_simul(order=1, irf=40) Y C I N W R m i pi A mu;