% forwardEM1D(sig, ztop, sep)
% ztop is a column vector containing the z-coordinate of the *top* of each layer [m]
% sig is a column vector containing the true electrical conductivity of each layer [S/m]
% orient is a scalar value specifying the dipole orientation (0 = vertical; 1 = horizontal)
% sepan  is a scalar value specifying the coil separation [m]

function ObsData = forwardEM1D(sig, ztop, orient, sepa)
    sigma_a = ones(length(sepa), 1);
    for i = 1:length(sepa)
        W = weightEM1D(ztop, sepa(i), orient(i)) ;
        sigma_a(i) = W*sig(:);
    end
    ObsData = [sigma_a sepa orient];
end