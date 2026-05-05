% RBC Model with Capital Adjustment Costs (Hayashi Framework)
% Author: Francisco Arizola Blua
var C I K N Y W R r q A;
varexo eps;

parameters beta delta alpha theta chi rho sigma phi;

% Calibration
beta  = 0.96;
delta = 0.025;
alpha = 0.36;
theta = 1;
chi   = 1;
rho   = 0.95;
sigma = 0.01;

% --- MACRO-PROCESADOR DE DYNARE ---
% Recibe el valor de phi desde MATLAB (0, 2 o 4).
@#ifndef phi_val
    @#define phi_val = 2
@#endif
phi = @{phi_val};
% ----------------------------------

% --------------------------------------------------------
% Steady State computation (Exact analytical solution)
% --------------------------------------------------------
A_ss = 1;
r_ss = 1/beta - 1;
R_ss = r_ss + delta;
q_ss = 1;

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

model;
% 1. Consumption Euler equation (Defines real interest rate r)
1/C = beta * (1/C(+1)) * (1 + r);

% 2. Tobin's q
q = 1 / (1 - phi*(I/K(-1) - delta));

% 3. q-Euler equation (Capital pricing)
% Nota timing: K_t+1 en math es el capital elegido hoy para mañana, o sea 'K' en Dynare.
q = beta * (C/C(+1)) * ( R(+1) + q(+1)*( 1 - delta + phi*(I(+1)/K - delta)*(I(+1)/K) - (phi/2)*(I(+1)/K - delta)^2 ) );

% 4. Labor supply
W/C = theta * N^chi;

% 5. Capital accumulation
K = (1 - delta)*K(-1) + I - (phi/2)*(I/K(-1) - delta)^2 * K(-1);

% 6. Goods market
Y = C + I;

% 7. Production function
Y = A * K(-1)^alpha * N^(1 - alpha);

% 8. Labor demand
W = (1 - alpha) * Y / N;

% 9. Capital demand
R = alpha * Y / K(-1);

% 10. TFP shock process
log(A) = rho * log(A(-1)) + eps;
end;

initval;
A = A_ss;
r = r_ss;
R = R_ss;
q = q_ss;
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

stoch_simul(order=1, irf=40) Y C N I W r q K;