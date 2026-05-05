% Baseline RBC Model
% Author: Francisco Arizola Blua
var C I K L Y W R A;
varexo eps;

parameters beta delta alpha gamma psi rho sigma;

% Calibration
beta  = 0.96;        % Discount factor
delta = 0.025;       % Depreciation rate
alpha = 0.36;        % Capital share
gamma = 1;           % Weight on consumption utility
psi   = 2;           % Weight on leisure utility
rho   = 0.95;        % Persistence of productivity
sigma = 0.01;        % Std. dev. of technology shock

% Steady state auxiliary ratios (outside Dynare model)
R_ss = 1 / beta - 1 + delta;
A_ss = 1;
K_L_ratio = (alpha * A_ss / R_ss)^(1 / (1 - alpha));
W_ss = (1 - alpha) * A_ss * K_L_ratio^alpha;
Y_L_ratio = A_ss * K_L_ratio^alpha;
I_L_ratio = delta * K_L_ratio;
C_L_ratio = Y_L_ratio - I_L_ratio;

L_ss = ( (gamma / psi) * (1 / C_L_ratio) * W_ss ) / ( 1 + (gamma / psi) * (1 / C_L_ratio) * W_ss );
K_ss = K_L_ratio * L_ss;
I_ss = I_L_ratio * L_ss;
Y_ss = Y_L_ratio * L_ss;
C_ss = C_L_ratio * L_ss;

% Model equations
model;
% 1. Euler equation
(1/C) = beta * (1/C(+1)) * (R(+1) + 1 - delta);
% 2. Intra-temporal condition (labor-leisure)
(psi/(1 - L)) = (gamma/C) * W;
% 3. Capital accumulation
K = (1 - delta) * K(-1) + I;
% 4. Goods market
Y = C + I;
% 5. Production function
Y = A * K(-1)^alpha * L^(1 - alpha);
% 6. Labor demand
W = (1 - alpha) * Y / L;
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
L = L_ss;
W = W_ss;
end;

steady;
check;

shocks;
var eps; stderr sigma;
end;

stoch_simul(order=1,irf=40) A Y L C I W R;