% function that build G using WeightFunction
% ztop = top layer position
% coilspacing
% orientation of the coils

function m = inversionEM2D(ztop, data, lamb, Wm, tol, itmax)
    ndata = size(data,1);
    nlay = length(ztop);
    G = zeros(ndata, nlay);   % initialize W matrix
    for i=1:size(G, 1)        % populate the matrix row-by-row
        R = weightEM2D(ztop, data(i,2), data(i,3));
        G(i, :) = R;
    end
    % setting up inverse model
    A = G'*G + lamb*Wm;
    b = G'*data(:, 1);
    % inversion
    m = cgs(A, b, tol, itmax);
end