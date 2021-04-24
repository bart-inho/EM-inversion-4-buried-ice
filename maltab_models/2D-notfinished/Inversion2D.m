% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; close all; clc; tic
tic
% load the dataset
load data2D

xlog = 0:0.2:10; %

% model parameters
ndata = size(data,1);
ztop = 0:0.1:5;
nlay = length(ztop);
sigma_a = data(:, 1);
nsteps = length(unique(data(:, 4)));

% noise
for i = 1:length(sigma_a)
    sigma_a(i) = sigma_a(i)+normrnd(0,1)*1e-5;
end
data(:, 1) = sigma_a;

ndata_p = size(data(data(:, 4) == 1), 1);
ndata = size(data,1);
nlay = length(ztop);

G_stack = cell(1, nlay);

for j = 1:nsteps
    G_prov = zeros(ndata_p, nlay);   % initialize W matrix
    for i=1:size(G_prov, 1)          % populate the matrix row-by-row
        R = weightEM2D(ztop, data(i,2), data(i,3));
        G_prov(i, :) = R;
    end
    G_stack{j} = G_prov;
end

G = blkdiag(G_stack{:});

% inversion

yolo = G';

A = G'*G;
b = G'*data(:, 1);

m = cgs(A, b, 1e-10, 1e2);
m = reshape(m, [], nsteps);

toc
figure()
inv = pcolor(xlog, ztop, m);
ylabel('depth [m]')
xlabel('width [m]')
set(gca, 'YDir','reverse')
shading interp;
saveas(inv, 'inversion_2d', 'png')