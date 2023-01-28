%% Run Tasks
% Runs and validates tasks

% Initial Orbit Frame Orientation Angles
% Hill/Body Frame: {i_r, i_theta, i_h}
% Circular orbit: w_b_n = theta_dot * i_h
PARAMS.w_b_n_lmo = [0; 0; 0.000884797]; % rad/s
PARAMS.w_b_n_gmo = [0; 0;  0.0000709003]; % rad/s
% Position vector on a circular orbit: r = radius* i_r
PARAMS.r_lmo = 3396.19 + 400; % km
PARAMS.r_gmo = 20424.2; % km
% (3-1-3) Euler Angles 
PARAMS.eu_lmo_init = [20;30;60];
PARAMS.eu_gmo_init = [0; 0; 250];





%% Task 1: Orbit Simulation

[n_pos_lmo, n_vel_lmo, ~] = simulate_orbit(PARAMS.r_lmo, PARAMS.eu_lmo_init, PARAMS.w_b_n_lmo, 1150, 1);
disp('rLMO and vLMO at 450s: ');
disp(n_pos_lmo(:,end))
disp(n_vel_lmo(:,end))
% Plotting
figure;
plot3(n_pos_lmo(1,:), n_pos_lmo(2,:), n_pos_lmo(3,:), 'o')
xlabel('n_1 (km)');
ylabel('n_2 (km)');
zlabel('n_3 (km)');
xl = xlim();
yl = ylim();
zl = zlim();
hold on;
line([0,2000], [0,0], [0,0], 'LineWidth', 3, 'Color', 'k');
line([0,0], [0,2000], [0,0], 'LineWidth', 3, 'Color', 'k');
line([0,0], [0,0], [0,2000], 'LineWidth', 3, 'Color', 'k');
title('Orbit Simulations LMO and GMO')

% Simulate GMO
[n_pos_gmo, n_vel_gmo, ~] = simulate_orbit(PARAMS.r_gmo, PARAMS.eu_gmo_init, PARAMS.w_b_n_gmo, 500000, 5000);
disp('rGMO and vGMO at 1150s: ');
disp(n_pos_gmo(:,end))
disp(n_vel_gmo(:,end))



% Plotting
hold on;
plot3(n_pos_gmo(1,:), n_pos_gmo(2,:), n_pos_gmo(3,:), 'o')

%% Task 2: Orbit Frame Orientation

% Find LMO DCM [HN] at 300s
[~, ~, dcm_n_b_t] = simulate_orbit(PARAMS.r_lmo, PARAMS.eu_lmo_init, PARAMS.w_b_n_lmo, 300, 1);

HN_300s = dcm_n_b_t(:,end-2:end);
disp('LMO HN Matrix at 300 s: ');
disp(HN_300s)

%% Task 3: Sun Pointing Reference Frame Orientation

dcm_rs_n = getRsN();
disp("Sun Pointing Reference Frame Orientation DCM: ")
disp(dcm_rs_n)
n_w_rs_n = [0; 0; 0];
% n_w_rs_n is 0, Rs frame moves constant with respect to the inertial frame

%% Task 4: Nadir-Pointing Reference Frame

% Return [RnN] and n_w_rn_n at t = 330s

[dcm_rn_n, n_w_rn_n] = getRnN(PARAMS.eu_lmo_init, PARAMS.w_b_n_lmo, 330, 1);


%% Task 5: GMO-Pointing Reference Frame

[dcm_rc_n] = getRcN(PARAMS, 330, 1);

n_w_rc_n = get_n_w_rc_n(PARAMS, 330, 0.01);


%% Task 6: Att Error Evaluation

% S/C initial conditions at t0
o_b_n = [0.3; -0.4; 0.5]; % S/C Initial attitude
b_w_b_n = [1.00; 1.75; -2.20]; % S/C Iniital body angular velocity

% Test calcAttError for sun-pointing at t0
[o_b_rs_n, b_w_b_rs_n] = calcAttErr(o_b_n, b_w_b_n, dcm_rs_n, n_w_rs_n);

% Test calcAttError for nadir-pointing at t0
[o_b_rn_n, b_w_b_rn_n] = calcAttErr(o_b_n, b_w_b_n, dcm_rn_n, n_w_rn_n);

% Test calcAttError for GMO reference pointing at t0
[o_b_rc_n, b_w_b_rc_n] = calcAttErr(o_b_n, b_w_b_n, dcm_rc_n, n_w_rc_n);


%% Task 7: Numerical Attitude Simulator

% Integrate State vector X forward 500s with u = 0
% Provide H = [I]*w_b_n at 500s expressed in the B frame

[tout, xout] = rk4(@sc_dynamics, b_w_b_n*pi/180, [0;0;0], 1, 500, 0);
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]

b_H = b_I * xout(:, end);

% Provide Rotational Kinetic Energy T = 1/2 w_b_n' * I * w_b_n at 500
% seconds
T = 1/2* xout(:,end)' * b_I * xout(:,end);

% Provide MRP attitude o_b_n (500s)
x0 = [o_b_n; b_w_b_n*pi/180];
[tout, xout] = rk4(@sc_dynamics_full, x0, [0;0;0], 1, 500, 1);
disp("MRP Attitude at 500s:")
disp(xout(1:3, end))
plot(tout, xout)

% Provide angular momentum vector n_H(500s) in inertial frame components
% Get DCM [BN] from the MRP
BN = mrp2dcm(xout(1:3,end));
n_H = BN' * b_H;
% If you apply a fixed control torque b_u = (0.01, -0.01, 0.02) Nm provide
% the attitude o_b_n (t = 100s)
[tout, xout] = rk4(@sc_dynamics_full, x0, [0.01;-0.01;0.02], 1, 100, 1);
disp("MRP Attitude with Fixed COntrol Torque at 100s:")
disp("b_u = [0.01; -0.01; 0.02]")
disp(xout(1:3,end))
plot(tout,xout)

%% Task 8: Sun Pointing Control

% Develop a full mission scenario function with mode switching
% First mode to develop is Sun Pointing
% S/C should go into sun pointing model immediately

[t, yout, tlm] = run_mission(400, 1, 'SUN-POINTING', PARAMS, @sc_dynamics_full);
figure;
plot(t, tlm.ctrl_err_att(1:3, :))
title('MRP Attitude Error BR')
grid on;
figure;
plot(t, yout(1:3, :))
title('MRP Attitude BN')
grid on;
figure;
plot(t, yout(4:6, :))
title('Angular Rate BN')
grid on;
figure;
plot(t, tlm.ctrl_err_rate(1:3, :))
title('Angular Rate Error')
grid on;

%% Task 9: Nadir Pointing Control

[t, yout, tlm] = run_mission(400, 1, 'NADIR-POINTING', PARAMS, @sc_dynamics_full);
figure;
plot(t, tlm.ctrl_err_att(1:3, :))
title('MRP Attitude Error BR')
grid on;
figure;
plot(t, yout(1:3, :))
title('MRP Attitude BN')
grid on;
figure;
plot(t, yout(4:6, :))
title('Angular Rate BN')
grid on;
figure;
plot(t, tlm.ctrl_err_rate(1:3, :))
title('Angular Rate Error')
grid on;

%% Task 10: GMO Pointing Control

[t, yout, tlm] = run_mission(400, 1, 'GMO-POINTING', PARAMS, @sc_dynamics_full);
figure;
plot(t, tlm.ctrl_err_att(1:3, :))
title('MRP Attitude Error BR')
grid on;
figure;
plot(t, yout(1:3, :))
title('MRP Attitude BN')
grid on;
figure;
plot(t, yout(4:6, :))
title('Angular Rate BN')
grid on;
figure;
plot(t, tlm.ctrl_err_rate(1:3, :))
title('Angular Rate Error')
grid on;


%% Task 11: Mission Scenario Simulation

% Propagate 6500 seconds 

% Sun is always in the n2 direction
% Control law must detumble satellite and either point the antenna (-b1)
% towards the GMO mother satellite, its sensor (b1) at Mars (nadir
% direction or -r direction) or its solar panels (b3) at the sun

% Due to high power demands of the telescope, the S/C needs to point its
% solar panels at the sun whenever the S/C is on the sunlit side of Mars
% (positive n2 position coordinate). In other words, when pointing at the
% sun, the S/C must point its solar panels axis b3 in the n2 direction

% To complete this 3D frame definition, assume that r1 must point in the
% -n1 direction

% When the S/C is on the shaded side (over the hemisphere with -n2 position
% coordinates) it must either enter communication or science mode. In
% science mode, the science platform axis b1 must point at the center of
% mars, thus in the Nadir direction. To complete the reference frame,
% assume that r2 must line up with the orbit along-track axis itheta

% The communication mode requires that the LMO and GMO S/C position vectors
% have an angular difference of less than 35 degrees. In this mode, the
% nano-satellite communication platform axis (-b1) must point towards the
% GMO S/C

% SC on sunlit side of mars: Point to Sun (b3 to sun)
% SC not on sunlit Mars side & GMO Visible: Point Antenna axis -b1 at GMO
% SC not on sunlt mars side & GMO is not visible: Point sensor axis b1
% along mars nadir direction