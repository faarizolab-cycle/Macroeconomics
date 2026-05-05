% RBC Model with Investment-Specific Technology (IST) Shocks
% Author: Francisco Arizola Blua
var C I K N Y W R A Z;
varexo eps_a eps_z;

parameters beta delta alpha theta chi rho_a rho_z sigma_a sigma_z;

% Calibration
beta    = 0.96;
delta   = 0.025;
alpha   = 0.36;
theta   = 1;
chi     = 1;       % Inverse Frisch elasticity
rho_a   = 0.95;
rho_z   = 0.95;    % Persistence of IST shock
sigma_a = 0.01;
sigma_z = 0.01;    % Std. dev. of IST shock

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
Z_ss = 1;
R_ss = 1/beta - 1 + delta;

% Ratios
K_N_ratio = (alpha * A_ss / R_ss)^(1 / (1 - alpha));
Y_N_ratio = A_ss * K_N_ratio^alpha;
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

model;
% 1. Euler equation
1/C = beta * (1/C(+1)) * (Z/Z(+1)) * (Z(+1)*R(+1) + 1 - delta);

% 2. Labor supply
theta * N^chi = W/C;

% 3. Capital accumulation
K = Z * I + (1 - delta) * K(-1);

% 4. Goods market (Resource constraint)
Y = C + I;

% 5. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% 6. Labor demand
W = (1 - alpha) * Y / N;

% 7. Capital demand
R = alpha * Y / K(-1);

% 8. TFP shock process
log(A) = rho_a * log(A(-1)) + eps_a;

% 9. IST shock process
log(Z) = rho_z * log(Z(-1)) + eps_z;
end;

initval;
A = A_ss;
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
var eps_a; stderr sigma_a;
var eps_z; stderr sigma_z;
end;

% Simulamos extrayendo las trayectorias para graficar
stoch_simul(order=1, irf=40) Y N C I W R A Z;