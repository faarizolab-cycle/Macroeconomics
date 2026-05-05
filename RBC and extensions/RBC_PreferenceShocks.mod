% RBC Model with Preference Shocks
% Author: Francisco Arizola Blua

var C I K N Y W R A PSI NU;
varexo eps_a eps_psi eps_nu;

parameters beta delta alpha theta chi rho_a rho_psi rho_nu sigma_a sigma_psi sigma_nu;

% Calibration
beta    = 0.96;
delta   = 0.025;
alpha   = 0.36;
theta   = 1;
chi     = 1;
rho_a   = 0.95;
rho_psi = 0.90;
rho_nu  = 0.90;
sigma_a   = 0.01;
sigma_psi = 0.01;
sigma_nu  = 0.01;

% Steady State (PSI=NU=1 in SS)
R_ss   = 1/beta - 1 + delta;
A_ss   = 1;
K_N    = (alpha * A_ss / R_ss)^(1/(1-alpha));
W_ss   = (1 - alpha) * A_ss * K_N^alpha;
Y_N    = A_ss * K_N^alpha;
I_N    = delta * K_N;
C_N    = Y_N - I_N;
N_ss   = (W_ss / (theta * C_N))^(1/(1+chi));
K_ss   = K_N * N_ss;
I_ss   = I_N * N_ss;
Y_ss   = Y_N * N_ss;
C_ss   = C_N * N_ss;

model;
% Euler equation (intertemporal shock PSI)
PSI / C = beta * (PSI(+1) / C(+1)) * (R(+1) + 1 - delta);

% Labor supply (intratemporal shock NU)
NU * theta * N^chi = W / C;

% Capital accumulation
K = (1 - delta) * K(-1) + I;

% Goods market
Y = C + I;

% Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% Labor demand
W = (1 - alpha) * Y / N;

% Capital demand
R = alpha * Y / K(-1);

% TFP process
log(A) = rho_a * log(A(-1)) + eps_a;

% Intertemporal preference shock
log(PSI) = rho_psi * log(PSI(-1)) + eps_psi;

% Intratemporal preference shock
log(NU) = rho_nu * log(NU(-1)) + eps_nu;
end;

initval;
A   = 1;
PSI = 1;
NU  = 1;
R   = R_ss;
K   = K_ss;
I   = I_ss;
Y   = Y_ss;
C   = C_ss;
N   = N_ss;
W   = W_ss;
end;

steady;
check;

shocks;
var eps_a;   stderr sigma_a;
var eps_psi; stderr sigma_psi;
var eps_nu;  stderr sigma_nu;
end;

stoch_simul(order=1, irf=40);