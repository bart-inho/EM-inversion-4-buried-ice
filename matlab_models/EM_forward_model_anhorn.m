clear all, close all, clc, 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////
tic

% initial conductivity
sig1 = 1e-3;
sig2 = 4e-3;
sig3 = 3e-3;

% frequency and depth
freq = 4e5;
depth = round(500/sqrt(sig1*freq)); % depth [m]

% building inital conductivity model
distance = 50; % A-B distance [m]
resmap = ones(depth, distance);
resmap(1:end,1:end) = sig1;
resmap(round(depth*(1/3)):round(depth*(2/3)),1:end) = sig2;
resmap(round(depth*(2/3)):end,1:end) = sig3;

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
% figure(2)
% hold on
% plot(resEM, centroid)
% ylabel('centroids')
% xlabel('conductivity \sigma')
% title('forward model')
% hold off

% /////////////////////////////////////////
% giving sigma_a
% ///////////////////////////////////////

coilspace = [3 4 5 6 7 8 9 10 12 15 20]; % coil spacing [m]
z1 = 10; z2 = 20; % depth [m]

sigma_a_h = sig1*(1-rh(z1, coilspace)) + sig2*(rh(z1, coilspace) - rh(z2, coilspace)) + sig3*rh(z2, coilspace);
% disp('horizontal sigma a = '); disp(sigma_a_h);
sigma_a_v = sig1*(1-rv(z1, coilspace)) + sig2*(rv(z1, coilspace) - rv(z2, coilspace)) + sig3*rv(z2, coilspace);
% disp('vertical sigma a = '); disp(sigma_a_v);
sigma_a_tot = [sigma_a_h sigma_a_v];
% disp('total sigma a = '); disp(sigma_a_tot)

% //////////////////////////////////////////////
% Inversion
% //////////////////////////////////////////////

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
legend('m0', 'homemade inversion', 'location', 'northwest')
xlabel('conductivity \sigma')
ylabel('depth [m]')
title('Homemade inversion')
set(gca, 'YDir','reverse')
hold off
saveas(inversion, 'inversion_fig.png')

% /////////////////////////////////
% 2D array combining 1D array
% /////////////////////////////////

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
% test using toolbox
% /////////////////////////////////

x = centroid';
y = sig3*(x<20) + sig2*(x>=20 & x < 30) + sig1*(x >= 30);

for i = 1:length(y)
    y(i) = y(i) + randn(1)*5e-4  ; % noise
end

model_SVM = fitrsvm(x,y, 'KernelFunction','gaussian','KernelScale','auto', 'Standardize',true);

sigma0 = 0.2;
kparams0 = [3.5, 6.2];
model_GMR = fitrgp(x,y,'KernelFunction','squaredexponential',...
     'KernelParameters',kparams0,'Sigma',sigma0);

yfit_GMR = predict(model_GMR,x);
yfit_SVM = predict(model_SVM,x);

% plot
% test = figure(4);
% hold on
% plot(y, x)
% plot(yfit_SVM, x)
% plot(yfit_GMR, x)
% % xlim([-0.001 0.015])
% % ylim([-5 55])
% legend('m0','SVM', 'GMR')
% xlabel('conductivity \sigma')
% ylabel('centroids')
% title('Using Machine Learning ToolBox')
% hold off
% saveas(test, 'test_fig.png')

% /////////////////////////////////
% functions 
% /////////////////////////////////

function Rh = rh(prof, coilspacing)
    Rh = 1./((4.*(prof./coilspacing).^2+1.^(1/2)-2.*(prof./coilspacing)));
end
function Rv = rv(prof, coilspacing)
    Rv = 1./((4.*(prof./coilspacing).^2+1.^(1/2)));
end
