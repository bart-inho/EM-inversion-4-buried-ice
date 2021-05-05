% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; clc; close all; tic
% load the dataset
load data2D

xlog = unique(data(:,4)); % horizontal size of the model
ztop = 0:0.05:5; % vertical size of the model
nx = length(xlog); % number of discretization layers horizontal
nz = length(ztop); % number of discretization layers vertical
ndata = size(data,1); % number of sigma_a

% noise
% Gaussian noise with standard deviation of X% of average data value
sigma_a = data(:, 1);
noi = 0.01;
sigma_a = sigma_a + noi*mean(abs(sigma_a))*randn(size(sigma_a));
data(:, 1) = sigma_a;

% setting up G matrix
ndata_p = size(data(data(:, 4) == 1), 1); % select data for each measurment point
G_stack = cell(1, nz); % predefine the list that will contain all 1D G marix

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

% set parameters for the regularization
% JI: alphax should generally be much bigger than alphaz for layered media
alphax = 1; % weight on model smoothness in x-direction
alphaz = 0.1; % weight on model smoothness in z-direction

% calculate data and model weighting matrices
[Dx,Dz] = smoothweightEM2D(nx,nz);
Wm = alphax*(Dx'*Dx) + alphaz*(Dz'*Dz);

% inversion
lamb = 1e-1; % noise fitting parameter
tol = 1e-10; % tolerance max
itmax = 1e4; % max iteration

A = G'*G + Wm*lamb; 
b = G'*data(:, 1);
[m1,fl1,rr1,it1,rv1] = cgs(A, b, tol, itmax);
m = reshape(m1, [], nx);

figure()
semilogy(0:length(rv1)-1,rv1/norm(b),'-')
hold on
yline(tol,'r--');
xlabel('Iteration number')
ylabel('Relative residual')

figure()
inv = pcolor(xlog, ztop, m);
title('Inverted model')
subtitle({['\lambda = ' num2str(lamb) ', \alpha_x = ' num2str(alphax), ', \alpha_z = ' num2str(alphaz)]...
    ['nx = ' num2str(nx) ', nz = ' num2str(nz) ', p_G(z) = ' num2str(noi)]...
    ['coilspacing = [ ' num2str(data(1, 2)) ', ' num2str(data(3, 2)) ', ' num2str(data(5, 2)) ' ]']})
set(gca, 'YDir','reverse')
xlabel('width [m]')
ylabel('depth [m]')
shading interp;
axis image
c = colorbar;
c.Label.String = '\sigma [S/m]';
% caxis([0 max(m1)]);
saveas(inv, 'inversion_2d', 'png')
toc