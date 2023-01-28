function [t, yout, tlm] = run_mission(tend, dt, mode, PARAMS, fcn)
% Simulate a mission using a runge-kutta 4 integrator

% Setup initial conditions
o_b_n = [0.3; -0.4; 0.5]; % S/C Initial attitude
b_w_b_n = [1.00; 1.75; -2.20]; % [deg/s] S/C Iniital body angular velocity
x0 = [o_b_n; b_w_b_n*pi/180];

t = 0:dt:tend;
yout = zeros(length(x0), length(t));
yout(:, 1)= x0; % Setup Initial Conditions

% Setup output telemtetry
tlm = struct();
tlm.target_mrp = zeros(3, length(t));
tlm.target_rate = zeros(3, length(t));
tlm.ctrl_err_att = zeros(3, length(t));
tlm.ctrl_err_rate = zeros(3, length(t));
tlm.u = zeros(3, length(t));
% Start mission simulation
for i = 1:(length(t)-1)
    % Check what mode we are in
    switch mode
        case 'SUN-POINTING'
            % Get Sun-Pointing Reference Frame
            % From Inertial to Reference
            dcm_r_n = getRsN(); % From Inertial to Sun-Pointing
            % Get desired angular velocity
            n_w_r_n = [0; 0; 0]; % rad/s
        case 'NADIR-POINTING'
            % Get Nadir-Pointing Reference Frame
            % From Inertial to Reference
            % Also get Desired Angular Velocity in inertial frame
            [dcm_r_n, n_w_r_n] = getRnN(PARAMS.eu_lmo_init, PARAMS.w_b_n_lmo, i, dt);
        case 'GMO-POINTING'
            % Get GMO-Pointing Reference Frame
            % From Inertial to Reference
            % Get Desired Angular Velocity in inertial frame
            dcm_r_n = getRcN(PARAMS, t, dt);
            n_w_r_n = get_n_w_rc_n(PARAMS, t, dt);
    end
    % Save setpoint variables to telemetry
    % Convert DCM to MRP to save in telemetry
    tlm.target_mrp(1:3, i) = dcm2mrp(dcm_r_n);
    tlm.target_rate(1:3, i) = n_w_r_n;
    % Calculate Control Error
    [o_b_r_n, b_w_b_r_n] = calcAttErr(yout(1:3,i), yout(4:6, i), dcm_r_n, n_w_r_n);
    % Save control error to telemetry
    tlm.ctrl_err_att(1:3, i) = o_b_r_n;
    tlm.ctrl_err_rate(1:3, i) = b_w_b_r_n;
    
    
    % Calculate output from PID controller
    b_u = pointing_controller(o_b_r_n, b_w_b_r_n);
    % Save control output to tlm
    tlm.u(1:3, i) = b_u;
    
    % Propogate dynamics
    k1 = fcn(yout(:, i), b_u, i);
    k2 = fcn(yout(:, i) + k1 * dt/2, b_u, i + dt/2);
    k3 = fcn(yout(:, i) + k2 * dt/2, b_u, i + dt/2);
    k4 = fcn(yout(:, i) + k3 * dt,   b_u, i + dt);
    
    % Calculate y(i+dt)
    yout(:, i+1) = yout(:, i) + (k1 + 2*k2 + 2*k3 + k4)*dt/6;
    
    % check norm and see if we're approaching singularity
    norm_o_b_n = norm(yout(1:3, i+1));
    if (norm_o_b_n > 1)
        yout(1:3, i+1) = -(yout(1:3,i+1) ./(norm_o_b_n^2));          
    end
    

    
    
end



end