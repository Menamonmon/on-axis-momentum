function [mrp] = dcm2mrp(dcm)
% Convert DCM to MRP representation
% DCM2Quat -> quat2mrp

Q = dcm2quat(dcm);
% Scalar number as its first column
mrp_1 = Q(2) / (1 + Q(1));
mrp_2 = Q(3) / (1 + Q(1));
mrp_3 = Q(4) / (1 + Q(1));

mrp = [mrp_1; mrp_2; mrp_3];




end