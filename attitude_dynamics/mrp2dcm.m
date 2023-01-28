function dcm = mrp2dcm(mrp)
% Convert MRP to DCM


% Skew-Sym Matrix
mrp_t = [0 -mrp(3) mrp(2); mrp(3) 0 -mrp(1); -mrp(2) mrp(1) 0];

mrp2 = mrp'*mrp;

dcm = eye(3) + ( 8 * (mrp_t)^2 - (4*(1-mrp2) * mrp_t) ) / (1+mrp2)^2;


end