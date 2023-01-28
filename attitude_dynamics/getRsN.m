function dcm_rs_n = getRsN()
% Get DCM from inertial to RS frame
% Rs frame = {-n1, cross(n2, -n1), n2}
% Euler Angles Sequence (313): {180, 90, 0}

dcm_rs_n = [-1 0 0; 0 0 1; 0 1 0];


end