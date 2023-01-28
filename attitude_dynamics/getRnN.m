function [dcm_rn_n, n_w_rn_n]  = getRnN(initial_angles, w_b_n, t, dt)
% Get From Inertial to Nadir-Pointing Reference Frame
% Gives RnN at given time t
% H = {ir, itheta, ih} Orbit Frame
% Rn = {-ir, itheta, cross(-ir, itheta)
%     - r1 points to planet
%     - r2 points to vel (RAM) direction

% From Orbit Frame to Nadir-Pointing Reference Frame 
% RnH = {0, 180, 0} euler angle 313 sequence
RnH = [-1 0 0; 0 1 0; 0 0 -1];


% Propogate orbit and get the dcm at time t\
% Get DCM From Inertial to Orbit Frame
[~, ~, dcm_b_n_t] = simulate_orbit(0, initial_angles, w_b_n, t, dt);
dcm_b_n_t = dcm_b_n_t(:,end-2:end); % Just need to get the last one
dcm_n_b_t = dcm_b_n_t';

% From inertial to Nadir pointing reference frame
% Compute [RnN] = [RnH] * [HN]
dcm_rn_n = RnH * dcm_b_n_t;

% Return angular velocity vector n_w_rn_n
n_w_rn_n = dcm_n_b_t * w_b_n;


end