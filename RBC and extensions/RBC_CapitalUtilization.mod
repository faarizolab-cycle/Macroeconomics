% RBC Model with Capital Utilization
% Author: Francisco Arizola Blua
var C I K N Y W R A u delta_u;
varexo eps;

parameters beta delta_0 alpha theta rho sigma phi_1 phi_2;

% Calibration Baseline
beta    = 0.96;
delta_0 = 0.025;
alpha   = 0.36;
theta   = 1;
rho     = 0.95;
sigma   = 0.01;

% --- MACRO-PROCESADOR DE DYNARE ---
% Recibe el valor de phi_2 desde MATLAB. Si no recibe nada, usa 0.1 por defecto.
@#ifndef phi2_val
    @#define phi2_val = 0.1
@#endif
phi_2 = @{phi2_val};
% ----------------------------------

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
u_ss = 1;
delta_ss = delta_0;

% Tasa de interés de estado estacionario
R_ss = 1/beta - 1 + delta_0;

% Restricción de calibración para que u_ss = 1
phi_1 = R_ss;

% Ratios
K_N_ratio = (alpha / R_ss)^(1 / (1 - alpha));
Y_N_ratio = K_N_ratio^alpha;
W_ss      = (1 - alpha) * Y_N_ratio;
I_N_ratio = delta_0 * K_N_ratio;

% Horas de trabajo en Estado Estacionario
N_ss = (1 - alpha) / ((1 - alpha) + theta * (1 - delta_0 * alpha / R_ss));

% Niveles
K_ss = K_N_ratio * N_ss;
Y_ss = Y_N_ratio * N_ss;
I_ss = I_N_ratio * N_ss;
C_ss = Y_ss - I_ss;

model;
% 1. Euler equation
1/C = beta * (1/C(+1)) * (R(+1)*u(+1) + 1 - delta_u(+1));

% 2. Labor supply
W/C = theta / (1 - N);

% 3. Optimal utilization
R = phi_1 + phi_2 * (u - 1);

% 4. Production function
Y = A * (u * K(-1))^alpha * N^(1 - alpha);

% 5. Labor demand
W = (1 - alpha) * Y / N;

% 6. Capital demand
R = alpha * Y / (u * K(-1));

% 7. Capital accumulation
K = (1 - delta_u) * K(-1) + I;

% 8. Goods market
Y = C + I;

% 9. Endogenous Depreciation
delta_u = delta_0 + phi_1 * (u - 1) + (phi_2 / 2) * (u - 1)^2;

% 10. TFP shock process
log(A) = rho * log(A(-1)) + eps;
end;

initval;
A = A_ss;
u = u_ss;
delta_u = delta_ss;
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
var eps; stderr sigma;
end;

stoch_simul(order=1, irf=40) Y N u C I W R;