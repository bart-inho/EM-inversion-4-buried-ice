% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles
tic
clear; close all; clc; tic

% load the dataset
load data1D
load TrueModel

% model parameters
ndata = size(data,1);
ztop = 0:0.1:2;
nlay = length(ztop);
sigma_a = data(:, 1);

% noise
for i = 1:length(sigma_a)
    sigma_a(i) = sigma_a(i)+normrnd(0,1)*5e-4;
end
data(:, 1) = sigma_a;

% smoothness
e = ones(nlay,1);
D = spdiags([e -2*e e], -1:1, nlay, nlay);

% regularization
Wm = D'*D;
lamb = 1e-2;

% inversion
m = inversionEM1D(ztop, data, lamb, Wm, 1e-10, 2e2);

% true model we know
sigma = ones(length(ztop), 1)*sig(1);
sigma(ceil(1/3*length(ztop)):ceil(2/3*length(ztop))) = sig(2);
sigma(ceil(2/3*length(ztop)):end) = sig(3);
toc

% plot
invers = figure();
plot(m,ztop)
hold on
plot(sigma, ztop)
legend('inverted model', 'initial model', 'location', 'best')
set(gca,'Ydir','reverse')
saveas(invers, 'inversion.png')