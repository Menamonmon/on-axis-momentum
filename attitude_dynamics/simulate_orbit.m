function [n_pos,n_vel, dcm_b_n_t] = simulate_orbit(radius, initial_angles, w_b_n, t, dt)

% Initialize variables

% Position vector on a circular orbit: r = radius* i_r
b_pos = [radius; 0; 0];
n_pos = [];
n_vel = [];
eu_angles = [];
eu_angles = [eu_angles initial_angles];
dcm_b_n_t = []; % struct to hold all of HN(t)

i = 1;
for i_t = 0:dt:t
    % Setup DCM (From Inertial to Body)
    dcm_b_n = get313DCM(eu_angles(:,i));
    dcm_b_n_t = [dcm_b_n_t dcm_b_n];
    
    % Update position vector
    % use DCM from body to inertial
    dcm_n_b = dcm_b_n';
    n_pos_val = dcm_n_b * b_pos;
    n_pos = [n_pos n_pos_val];
    
    % Update velocity vector
    b_vel_val = cross(w_b_n, b_pos);
    n_vel_val = dcm_n_b * b_vel_val;
    n_vel = [n_vel n_vel_val];

    
    
    % Get new euler rates (From body to inertial)
    
    % From Inertial to Body
    dcm_rate_b_n = get313DiffEq(eu_angles(:,i));
    % Get From Body to Inertial
    dcm_rate_n_b = dcm_rate_b_n';
    eul_rates = dcm_rate_n_b * w_b_n; % From body to inertial
    % Theta = Theta_initial + t * theta_dot
    new_angles_rad = eu_angles(:,i) * (pi/180) + eul_rates * dt;
    

    eu_angles = [eu_angles new_angles_rad*(180/pi)];
    i = i + 1;
        
end



end

