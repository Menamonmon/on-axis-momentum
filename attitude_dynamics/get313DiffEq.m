function matrix = get313DiffEq(eu_angles)

% Gets Kinematic Differential Equations matrix for (3-1-3) Euler angle
% sequence. From inertial frame to body frame

theta_2 = eu_angles(2) * (pi/180);
theta_3 = eu_angles(3) * (pi/180);

matrix = [sin(theta_3)*sin(theta_2) cos(theta_3) 0; 
          cos(theta_3)*sin(theta_2) -sin(theta_3) 0; 
          cos(theta_2) 0 1];


end