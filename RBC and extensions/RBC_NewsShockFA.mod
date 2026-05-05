% RBC Model with News Shocks
% Author: Francisco Arizola Blua
var C I K N Y W R Z;
varexo eps_surprise eps_news;

parameters beta delta alpha psi rho_z sigma_surp sigma_news;

% Calibration
beta  = 0.96;
delta = 0.025;
alpha = 0.36;
psi   = 2;
rho_z = 0.95;
sigma_surp = 0.01;
sigma_news = 0.01;

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
Z_ss = 1;
R_ss = 1/beta - 1 + delta;

% Exact SS for Labor (N)
N_ss = ( (1 - alpha) / psi ) / ( 1 - delta * alpha / R_ss + (1 - alpha) / psi );

% Ratios and Levels
K_N_ratio = (alpha / R_ss)^(1 / (1 - alpha));
K_ss = K_N_ratio * N_ss;
Y_ss = K_N_ratio^alpha * N_ss;
I_ss = delta * K_ss;
C_ss = Y_ss - I_ss;
W_ss = (1 - alpha) * Y_ss / N_ss;

model;
% 1. Euler equation
1/C = beta * (1/C(+1)) * (R(+1) + 1 - delta);

% 2. Labor supply
(psi * C) / (1 - N) = W;

% 3. Capital accumulation
K = (1 - delta) * K(-1) + I;

% 4. Goods market
Y = C + I;

% 5. Production function
Y = Z * K(-1)^alpha * N^(1 - alpha);

% 6. Labor demand
W = (1 - alpha) * Y / N;

% 7. Capital demand
R = alpha * Y / K(-1);

% 8. TFP shock process
% eps_surprise: hits today
% eps_news(-8): was announced 8 periods ago, hits today. 
% Therefore, an eps_news shock today will hit Z in 8 periods.
log(Z) = rho_z * log(Z(-1)) + eps_surprise + eps_news(-8);
end;

initval;
Z = Z_ss;
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
var eps_surprise; stderr sigma_surp;
var eps_news;     stderr sigma_news;
end;

% Simulamos un horizonte de 40 periodos para ver bien el efecto de anticipación (8 periodos)
stoch_simul(order=1, irf=40) Z Y N C I W R;