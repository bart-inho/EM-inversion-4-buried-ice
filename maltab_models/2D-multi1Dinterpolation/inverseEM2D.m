% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; close all; clc; tic

% load the dataset
load data2D

xlog = 0:1:10; %

% model parameters
ndata = size(data,1);
ztop = 0:0.1:2;
nlay = length(ztop);
sigma_a = data(:, 1);

figure()
plot(sigma_a, '*')
hold on

% smallness
I = eye(nlay);

% smoothnesse = ones(nlay,1);
e = ones(nlay,1);
D = spdiags([e -2*e e], -1:1, nlay, nlay);
D(1, :) = 0; D(end, :) = 0;

% regularization parameters
alpha_s = 1 ;
alpha_z = 1 ;

% regularization
Wm = alpha_s*I + alpha_z*(D'*D);
lamb = 1e-6;

% noise
for i = 1:length(sigma_a)
    sigma_a(i) = sigma_a(i)+normrnd(0,1)*1e-5;
end
data(:,1) = sigma_a;

plot(sigma_a, 'x')
hold off

m = zeros(length(ztop),length(xlog));

for i = xlog
    m(:, i+1) = Inversion2D(ztop, data(data(:, 4)==i,1:3), lamb, Wm, 1e-10, 2e2);
end

figure()
inv = pcolor(xlog, ztop, m);
ylabel('depth [m]')
xlabel('width [m]')
shading interp;
saveas(inv, 'inversion_2d', 'png')