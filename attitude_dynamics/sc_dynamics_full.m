function x_dot = sc_dynamics_full(x, u, t)
% x = [o_b_n b_w_b_n]
% Calculate x_dot at time t
b_I = [10 0 0; 0 5 0; 0 0 7.5]; % [kg-m^2]

% Get States
o_b_n = x(1:3); 
b_w_b_n = x(4:6);

% Angular acceleration differential equation
b_w_b_n_skew = [0          -b_w_b_n(3) b_w_b_n(2); 
                b_w_b_n(3)  0         -b_w_b_n(1); 
               -b_w_b_n(2)  b_w_b_n(1) 0];


b_w_dot_b_n = b_I \ ((-b_w_b_n_skew * (b_I * b_w_b_n)) + u);


% MRP Differential Equation
skew_o_b_n = [0        -o_b_n(3) o_b_n(2); 
              o_b_n(3)  0       -o_b_n(1); 
             -o_b_n(2)  o_b_n(1) 0];


o_dot_b_n = 1/4 * ( (1-(o_b_n.'*o_b_n)) * eye(3) + 2 * skew_o_b_n + 2*(o_b_n*o_b_n.')) * b_w_b_n;

% Output x_dot
x_dot = [o_dot_b_n; b_w_dot_b_n];

end