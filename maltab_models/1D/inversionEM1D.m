% Function to invert a dataset
% ztop : vertical discretization
% data : data from field or CreateData code
% lamb : regularization parameter
% Wm : Smoothness
% tol : min tolerange for cgs function
% itmax : max iteration to find minimum

function m = inversionEM1D(ztop, data, lamb, Wm, tol, itmax)
    ndata = size(data,1);
    nlay = length(ztop);
    G = zeros(ndata, nlay);   % initialize W matrix
    for i=1:size(G, 1)        % populate the matrix row-by-row
        R = weightEM1D(ztop, data(i,2), data(i,3));
        G(i, :) = R;
    end
    size(G)
    % setting up inverse model
    A = G'*G + lamb*Wm;
    b = G'*data(:, 1);
    % inversion
    m = cgs(A, b, tol, itmax);
end