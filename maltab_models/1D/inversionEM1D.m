% Function to invert a dataset
% ztop : vertical discretization
% data : data from field or CreateData code
% lamb : regularization parameter
% Wm : Smoothness
% tol : min tolerange for cgs function
% itmax : max iteration to find minimum

function [m, G, A] = inversionEM1D(nz, data, lamb, Wm, Wd, m0, tol, itmax)
    ndata = size(data,1);
    nlay = length(nz);
    G = zeros(ndata, nlay);   % initialize W matrix
    for i=1:size(G, 1)        % populate the matrix row-by-row
        R = weightEM1D(nz, data(i,2), data(i,3));
        G(i, :) = R;
    end
    % setting up inverse model
    A = G'*Wd*G + lamb*Wm;
    b = G'*Wd*data(:,1) + lamb*Wm*m0;
    % inversion
    m = lsqr(A, b, tol, itmax);
end