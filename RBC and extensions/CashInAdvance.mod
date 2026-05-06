% RBC Model with Cash-in-Advance (CIA) Constraint
% Author: Francisco Arizola Blua
var C I K N Y W R xi i pi m g A;
varexo eps_a eps_m;

parameters beta delta alpha theta rho_a rho_m g_bar sigma_a sigma_m;

% Calibration
beta    = 0.96;
delta   = 0.025;
alpha   = 0.36;
theta   = 1;
rho_a   = 0.95;
rho_m   = 0.50;      % Persistence of money growth shock
g_bar   = 1.01;      % Steady-state gross money growth rate (~1% inflation)
sigma_a = 0.01;
sigma_m = 0.005;

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss  = 1;
g_ss  = g_bar;
pi_ss = g_bar - 1;

R_ss = 1/beta - 1 + delta;
i_ss = g_bar/beta - 1;

% Ratios
K_N_ratio = (alpha * A_ss / R_ss)^(1 / (1 - alpha));
Y_N_ratio = A_ss * K_N_ratio^alpha;
W_ss      = (1 - alpha) * Y_N_ratio;
I_N_ratio = delta * K_N_ratio;
C_N_ratio = Y_N_ratio - I_N_ratio;

% Exact SS for Labor (N) with CIA distortion
N_ss = (beta * W_ss) / (theta * g_ss * C_N_ratio + beta * W_ss);

% Levels
K_ss = K_N_ratio * N_ss;
Y_ss = Y_N_ratio * N_ss;
I_ss = I_N_ratio * N_ss;
C_ss = C_N_ratio * N_ss;

% Multiplier and Money Steady States
xi_ss = (1/C_ss) * (1 - beta/g_ss);
m_ss  = C_ss * g_ss;

model;
% 1. Euler equation for capital
1/C - xi = beta * (1/C(+1) - xi(+1)) * (R(+1) + 1 - delta);

% 2. Euler equation for bonds (Fisher equation)
1/C - xi = beta * (1/C(+1) - xi(+1)) * (1 + i) / (1 + pi(+1));

% 3. FOC for money
1/C - xi = beta * (1/C(+1)) / (1 + pi(+1));

% 4. Labor supply (Notice the logarithmic preference over leisure)
theta / (1 - N) = (1/C - xi) * W;

% 5. Cash-in-Advance (CIA) Constraint (binding)
m(-1) / (1 + pi) = C;

% 6. Capital accumulation
K = (1 - delta) * K(-1) + I;

% 7. Resource constraint
Y = C + I;

% 8. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% 9. Labor demand
W = (1 - alpha) * Y / N;

% 10. Capital demand
R = alpha * Y / K(-1);

% 11. Money growth identity (in real terms)
g = (m / m(-1)) * (1 + pi);

% 12. TFP shock process
log(A) = rho_a * log(A(-1)) + eps_a;

% 13. Money growth process
log(g) = (1 - rho_m)*log(g_bar) + rho_m*log(g(-1)) + eps_m;
end;

initval;
A  = A_ss;
g  = g_ss;
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
xi = xi_ss;
end;

steady;
check;

shocks;
var eps_a; stderr sigma_a;
var eps_m; stderr sigma_m;
end;

stoch_simul(order=1, irf=40) Y C I N W R m i pi A g;