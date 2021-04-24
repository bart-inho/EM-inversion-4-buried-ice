% forwardEM1D(sig, ztop, sep)
% ztop is a column vector containing the z-coordinate of the *top* of each layer [m]
% sigma is a column vector containing the true electrical conductivity of each layer [S/m]
% orient is a scalar value specifying the dipole orientation (0 = vertical; 1 = horizontal)
% coilsep is a scalar value specifying the coil separation [m]

function ObsData = forwardEM2D(sig, ztop, sepa, orient, x)
    x = ones(length(sepa), 1)*x;
    sigma_a = ones(length(sepa), 1);
    for i = 1:length(sepa)
        W = weightEM2D(ztop, sepa(i), orient(i)) ;
        sigma_a(i) = W*sig(:);
    end
    ObsData = [sigma_a sepa orient x];
end