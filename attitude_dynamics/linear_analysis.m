% Linear Analysis


% Open Loop Linearized Equations of Motion
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]

% x = [o_b_n b_w_b_n]
% x = [o1 o2 o3 w1 w2 w3]

A = [0 0 0 1/4 0 0; 0 0 0 0 1/4 0; 0 0 0 0 0 1/4; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
B = [zeros(3); inv(b_I)];
C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
D = [zeros(3); zeros(3)];

sys_ol = ss(A, B, C, D);

step(sys_ol)
pole(sys_ol); % 6 poles all at the origin

% PD-Like Nonlinear Control Law

% Time for the tracking errors should be 1/e the original size) = 120
% seconds

% All decay time constants should be 120 seconds or less

% Closed loop resonse for attitude should be either crit damped or
% under-damped. At least be crit damped with eta = 1, while other modes
% will have eta <=1

% Solve for P and K
T = 120;
eta = 1;
P = (2 * b_I)/T;
K = 1/b_I(3,3) * (P(3,3) / eta)^2;
% Verify this K gain keeps <=1 damping requirement for the other modes
eta_1 = P(1,1) / sqrt(K * b_I(1,1));
eta_2 = P(2,2) / sqrt(K * b_I(2,2));
eta_3 = P(3,3) / sqrt(K * b_I(3,3));
A_cl = [0 0 0 1/4 0 0; 0 0 0 0 1/4 0; 0 0 0 0 0 1/4; -K*inv(b_I) -P*inv(b_I)]; 
B_cl = [zeros(3); zeros(3)];
C_cl = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
D_cl = [zeros(3); zeros(3)];

sys_cl = ss(A_cl, B_cl, C_cl, D_cl);
o_b_n = [0.3; -0.4; 0.5]; % S/C Initial attitude
b_w_b_n = [1.00; 1.75; -2.20]; % S/C Iniital body angular velocity
x0 = [o_b_n; b_w_b_n*pi/180];
% Response to initial conditions
[y, tout, x] = lsim(sys_cl, zeros(3,4001), 0:0.1:400, x0);

figure;
plot(tout, y(:, 1:3))
title('MRP BN')
figure;
plot(tout, y(:, 4:6))
title('Angular Rates')