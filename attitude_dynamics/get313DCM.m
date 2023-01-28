function dcm_matrix = get313DCM(eu_angles)
% Get 3-1-3 DCM from inertial to body 

% dcm = angle2dcm(rot1, rot2, rot3, 'ZXZ')
theta_1 = eu_angles(1) * (pi/180);
theta_2 = eu_angles(2) * (pi/180);
theta_3 = eu_angles(3) * (pi/180);

r_1 = [cos(theta_1) sin(theta_1) 0; -sin(theta_1) cos(theta_1) 0; 0 0 1];
r_2 = [1 0 0; 0 cos(theta_2) sin(theta_2); 0 -sin(theta_2) cos(theta_2)];
r_3 = [cos(theta_3) sin(theta_3) 0; -sin(theta_3) cos(theta_3) 0; 0 0 1];

dcm_matrix = r_3 * r_2 * r_1;


end