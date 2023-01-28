function [RcN] = getRcN(PARAMS,t, dt)
% Get GMO-Pointing Reference Frame Orientation
%
% Position LMO S/C is away from GMO S/C:
%   dr = r_lmo - r_gmo
%
% Rc = {-dr/|dr|, cross(dr, n3 / |cross(dr, nr)|, cross(r1,r2)}
%       - -r1 points to GMO S/C
%       - r2 fully defines a 3-D reference frame
%       - r3 fully defines a 3-D reference frame

r_lmo = PARAMS.r_lmo;
r_gmo = PARAMS.r_gmo;
w_b_n_lmo = PARAMS.w_b_n_lmo;
w_b_n_gmo = PARAMS.w_b_n_gmo;
eu_lmo_init = PARAMS.eu_lmo_init;
eu_gmo_init = PARAMS.eu_gmo_init;


% Get position of LMO S/C wrt to inertial frame at time t
[n_pos_lmo, ~, ~] = simulate_orbit(r_lmo, eu_lmo_init, w_b_n_lmo, t, dt);


% Get position of GMO S/C wrt to inertial frame at time t
[n_pos_gmo, ~, ~] = simulate_orbit(r_gmo, eu_gmo_init, w_b_n_gmo, t, dt);

dr = n_pos_gmo(:,end) - n_pos_lmo(:,end);
r1 = -(dr) / norm(dr);
r2 = cross(dr, [0 0 1]') / norm(cross(dr, [0 0 1]'));
r3 = cross(r1, r2);

RcN = [r1 r2 r3]';

end


