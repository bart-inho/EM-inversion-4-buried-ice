clear all, close all, clc, 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////
tic

% initial conductivity
sig1 = 1e-3;
sig2 = 4e-3;
sig3 = 3e-3;
sig = [sig1; sig2; sig3];

% frequency and depth
freq = 4e5;
depth = round(500/sqrt(sig1*freq)); % depth [m]
depth_layer = [round(depth*1/3); round(depth*2/3)];
% building inital conductivity model
distance = 50; % A-B distance [m]
resmap = ones(depth, distance);
resmap(1:end,1:end) = sig1;
resmap(depth_layer(1):depth_layer(2),1:end) = sig2;
resmap(depth_layer(2):end,1:end) = sig3;

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

coil = 1;
ori = 0;

sigma_a = forwardEM1D(depth_layer, sig, ori, coil);
disp(sigma_a)

% /////////////////////////////////////////////
% Inversion using sigma_a
% /////////////////////////////////////////////

% to finish

% //////////////////////////////////////////////
% Polynomial regression
% //////////////////////////////////////////////
% Invented datas, without taking sigma_a in count (in progress)
% setting up x and y
x = resEM';
y = centroid;
for i = 1:length(x)
    x(i) = x(i) + randn(1)*5e-5 ;
end % noise

tic
% setting up inversion functions
k = 0:1:8;
G = zeros(length(y), length(k));
G(:, :) = y(:).^repmat(k, length(y), 1);
beta = inv(G'*G)*(G'*x);
y1=zeros(length(y), 1);

% polynomial regression
for i = 1:length(k)
   y1(:) = y1(:) + beta(i)*y(:).^k(i);
end 
toc

% plot
inversion = figure(3);
hold on
plot(x, y)
plot(y1, y, 'r')
% xlim([-0.001 0.015])
% ylim([-5 55])
legend('m0', 'homemade regression', 'location', 'northwest')
xlabel('conductivity \sigma')
ylabel('depth [m]')
title('Homemade inversion')
set(gca, 'YDir','reverse')
hold off
saveas(inversion, 'inversion_fig.png')

% /////////////////////////////////
% 2D array combining 1D array
% /////////////////////////////////
% represent regression in 2D

measurement_spacing = thkness; % [m]
terrain_size = distance; % [m]
dist = ones(1, terrain_size)*measurement_spacing; % 

model2d = dist.*y1 ;

inversion_2d = figure(5);
pcolor(model2d)
xlabel('x [m]')
ylabel('depth [m]')
title('conductvity model inverted')
colorbar
saveas(inversion_2d, 'inversion_2d.png')

% /////////////////////////////////

toc

% /////////////////////////////////
% functions 
% /////////////////////////////////

% -> weight W could be better 

function C = forwardEM1D(ztop, sigma, orient, sep)
    if orient == 1 % 1 = horizontal
        R = ((4.*((ztop./sep).^2)+1).^(1/2))-2.*(ztop./sep);
    else           % 0 = vertical
        R = ((4.*((ztop./sep).^2)+1).^(1/2)).\1;
    end
    C = sigma(1)*(1-R(1)) + sigma(end)*R(end);
        for j = 2:length(sigma)-1
           C = C + sigma(j).*(R(j-1) - R(j));
        end 
end