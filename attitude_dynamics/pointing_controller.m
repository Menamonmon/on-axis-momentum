function b_u = pointing_controller(o_b_r, b_w_b_r)

% Calculate control torques u (x, y, and z) 
% Controller to drive body frame B towards reference R
K = 0.0021;
P = [0.1667 0 0;0 0.0833 0; 0 0 0.1250];
%K = 0.0056;
%P = 0.1667;
b_u = (-K * o_b_r) + (-P * b_w_b_r);


end