% z = top layer height
% sep = coil spacing
% ori = coil orientation (0, 1)
function Weight = Weight1D(z, sep, ori)
    if ori == 0
        w = (4*(z./sep).^2+1).^(-1/2);             % vertical dipole
    elseif ori == 1
        w = (4*(z./sep).^2+1).^(1/2) - 2*z./sep;    % horizontal dipole
    end
    Weight = w(1:end-1)-w(2:end); % R = weight (R(zn) - R(zn+1))
    Weight(end+1) = w(end); % R(end) = only R(z)
    Weight = Weight';
end