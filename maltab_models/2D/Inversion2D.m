% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; clc; tic
tic
% load the dataset
load data2D

xlog = unique(data(:,4)); % horizontal size of the model
ztop = 0:0.1:10; % vertical size of the model

nx = length(xlog); % number of discretization layers horizontal
nz = length(ztop); % number of discretization layers vertical

ndata = size(data,1); % number of sigma_a

% noise
% Gaussian noise with standard deviation of X% of average data value
sigma_a = data(:, 1);
sigma_a = sigma_a + 0.01*mean(abs(sigma_a))*randn(size(sigma_a));
data(:, 1) = sigma_a;

ndata_p = size(data(data(:, 4) == 1), 1); % select data for each measurment point

G_stack = cell(1, nz); % predefine the list that will contain all 1D G marix

for j = 1:nx
    G_prov = zeros(ndata_p, nz);   % initialize W matrix
    for i=1:size(G_prov, 1)          % populate the matrix row-by-row
        R = weightEM2D(ztop, data(i,2), data(i,3));
        G_prov(i, :) = R;
    end
    G_stack{j} = G_prov;
end

G = sparse(blkdiag(G_stack{:})); % stack, blockdiag and sparse

% set parameters for the regularization
% JI: alphax should generally be much bigger than alphaz for layered media
alphax = 1; % weight on model smoothness in x-direction
alphaz = 0.1; % weight on model smoothness in z-direction

% calculate data and model weighting matrices
[Dx,Dz] = smoothweightEM2D(nx,nz);
Wm = alphax*(Dx'*Dx) + alphaz*(Dz'*Dz);


% inversion
lamb = 1e-2;
A = G'*G + Wm*lamb; 
b = G'*data(:, 1);
m = cgs(A, b, 1e-10, 1000);
m = reshape(m, [], nx);

toc

figure()
inv = pcolor(xlog, ztop, m);
title(['\lambda = ' num2str(lamb) ', \alpha_x = ' num2str(alphax) ', \alpha_z = ' num2str(alphaz)]) 
ylabel('depth [m]')
xlabel('width [m]')
set(gca, 'YDir','reverse')
colorbar
% shading interp;
saveas(inv, 'inversion_2d', 'png')