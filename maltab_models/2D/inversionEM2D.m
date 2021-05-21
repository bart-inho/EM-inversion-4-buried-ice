function [m, m1, G, A] = inversionEM2D(ztop, nz, nx, data, lamb, Wm, Wd, m0, tol, itmax)
% setting up G matrix
ndata_p = size(data(data(:, 4) == 1), 1); % select data for each measurment point
G_stack = cell(1, nz); % predefine the list that will contain all 1D G marix
d = data(:,1);
for j = 1:nx
    G_prov = zeros(ndata_p, nz);   % initialize W matrix
    for i=1:size(G_prov, 1)        % populate the matrix row-by-row
        R = weightEM2D(ztop, data(i,2), data(i,3));
        G_prov(i, :) = R;
    end
    G_stack{j} = G_prov;
end

% G matrix
G = sparse(blkdiag(G_stack{:})); % stack, blockdiag and sparse

% solve linear system  %JI: note additions here to reflect above changes
A = G'*Wd*G + lamb*Wm; 
b = G'*Wd*d + lamb*Wm*m0;
m1 = cgs(A, b, tol, itmax);
m = reshape(m1,nz,nx);
