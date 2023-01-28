function x_dot = sc_dynamics(x, u, t)
% x = [w_b_n(1) w_b_n(2) w_b_n(3)]
% Calculate x_dot at time t
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]


b_w_b_n = x;

% Skew Matrix b_w_b_n_skew
b_w_b_n_skew = [0          -b_w_b_n(3) b_w_b_n(2); 
                b_w_b_n(3)  0         -b_w_b_n(1); 
               -b_w_b_n(2)  b_w_b_n(1) 0];


b_w_dot_b_n = inv(b_I)* (-b_w_b_n_skew * b_I * b_w_b_n + u);

x_dot = b_w_dot_b_n;



end