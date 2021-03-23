clear all, close all, clc, 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////
tic

% initial conductivity
sig = [4e-3; 2e-3; 8e-3];

% frequency and depth
freq = 4e5;
depth = round(500/sqrt(sig(1)*freq)); % depth [m]

z = [round(depth*1/3); round(depth*2/3)];

% building inital conductivity model
distance = 50; % A-B distance [m]
resmap = ones(depth, distance);
resmap(1:end,1:end) = sig(1);
resmap(z(1):z(2),1:end) = sig(2);
resmap(z(2):end,1:end) = sig(3);

% plotting the initial model 
initial_model = figure(1);
pcolor(resmap)
xlabel('x')
ylabel('depth')
title('conductvity model')
colorbar
saveas(initial_model, 'initial_model.png')

% one dimension model :
position = 22; % random x on the 2D model to go in the 1D model
nlay = depth; % 
thkness = 1; % [m]
thk = ones(1, nlay); % creating a thickness matrix
thk = thk*thkness; % thickness matrix
centroid = cumsum(thk) - thk/2; % layer centroids
resEM = ones(1, nlay)*resmap(1, position); % resistivity model

% indexing "inversion" (0:0 must be at the reference)
for i = centroid(2:end)
    if i < depth
        ab = depth+1-round(i);
        ab2 = floor(i/thkness);
        resEM(ab2) = resmap(ab, position);
    end
end

% plotting the forward model
figure(2)
hold on
plot(resEM, centroid)
ylabel('centroids')
xlabel('conductivity \sigma')
title('forward model')
hold off

% /////////////////////////////////////////
% giving sigma_a
% ///////////////////////////////////////

% forwardEM1D(ztop, sigma, orient, sep)
% ztop is a column vector containing the z-coordinate of the *top* of each layer [m]
% sigma is a column vector containing the true electrical conductivity of each layer [S/m]
% orient is a scalar value specifying the dipole orientation (0 = vertical; 1 = horizontal)
% sep is a scalar value specifying the coil separation [m]
% if two direction desired : sigma_a = [sigma_a_v sigma_a_h] ?

coil = 1;
ori = 0;
sigma_a = forwardEM1D(z, sig, ori, coil);
disp(sigma_a)

% /////////////////////////////////////////////
% Inversion using sigma_a
% /////////////////////////////////////////////

% ??

toc

% /////////////////////////////////
% functions 
% /////////////////////////////////

% -> weight W could be better 

function C = forwardEM1D(ztop, sigma, orient, sep)
    if orient == 1     % 1 = horizontal
        R = ((4.*((ztop./sep).^2)+1).^(1/2))-2.*(ztop./sep);
    elseif orient == 0 % 0 = vertical
        R = ((4.*((ztop./sep).^2)+1).^(1/2)).\1;
    end
    C = sigma(1)*(1-R(1)) + sigma(end)*R(end);
        for j = 2:length(sigma)-1
           C = C + sigma(j).*(R(j-1) - R(j));
        end
end