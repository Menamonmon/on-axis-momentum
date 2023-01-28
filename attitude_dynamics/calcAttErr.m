function [o_b_r, b_w_b_r] = calcAttErr(o_b_n, b_w_b_n, RN, n_w_r_n)
% Calculate Attitude and Angular velocity tracking error of current body
% frame B relative to reference frame R

BN = mrp2dcm(o_b_n);
NR = RN';
o_b_r = dcm2mrp(BN*NR);
b_w_b_r = (b_w_b_n - BN * n_w_r_n);




end