function n_w_rc_n = get_n_w_rc_n(PARAMS, t, dt)
% Calculate inertial angular velocity of Rc wrt inertial frame 


% Get RcN at time t
RcN_t = getRcN(PARAMS, t, dt);
% Get RcN at time t + dt
RcN_t_dt = getRcN(PARAMS, t+dt, dt);

skew_w = ((RcN_t_dt - RcN_t) / dt) * RcN_t';

% Deskew

n_w_rc_n = [-skew_w(3, 2); -skew_w(1, 3); -skew_w(2, 1)];



end