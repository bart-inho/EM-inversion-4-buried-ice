clear, close, clc, 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////

% building initial conductivity model
depth = 30; % depth
resmap = ones(30, 50);
sig1 = 1e-3;
sig2 = 1e-4;
sig3 = 2e-4;
resmap(1:end,1:end) = sig1;
resmap(10:20,1:end) = sig2;
resmap(20:end,1:end) = sig3;

% plotting the initial model 
figure(1)
pcolor(resmap)
xlabel('x')
ylabel('depth')
title('conductvity model')
colorbar

% one dimension model :
position = 22; % random x on the 2D model to go in the 1D model
nlay = 100; % n arbitrary layer
thkness = 0.5; % [m]
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

coilspace = [3 4 5 6 7 8 9 10 12 15 20]; % coil spacing [m]
z1 = 10; z2 = 20; % depth [m]

sigma_a_h = sig1*(1-rh(z1, coilspace)) + sig2*(rh(z1, coilspace) - rh(z2, coilspace)) + sig3*rh(z2, coilspace);
% disp('horizontal sigma a = '); disp(sigma_a_h);
sigma_a_v = sig1*(1-rv(z1, coilspace)) + sig2*(rv(z1, coilspace) - rv(z2, coilspace)) + sig3*rv(z2, coilspace);
% disp('vertical sigma a = '); disp(sigma_a_v);

sigma_a_tot = [sigma_a_h sigma_a_v];
disp('total sigma a = '); disp(sigma_a_tot)

% //////////////////////////////////////
% Model m0
% //////////////////////////////////////

% m0 = ones(1, length(resEM))*sigma_a_h;
 m0 = resEM;
for i = 1:length(m0)
    m0(i) = m0(i) + randn(1)*1.5e-4    ; % no noise
end

hold on
plot(m0, centroid, 'x')
legend('true model', 'm0 model')
% ylim([-1e-4 12e-4])
hold off

% //////////////////////////////////////////////
% Inversion
% ////////////////////////////////////////////




% /////////////////////////////////
% functions 
% /////////////////////////////////
function Rh = rh(prof, coilspacing)
    Rh = 1./((4.*(prof./coilspacing).^2+1.^(1/2)-2.*(prof./coilspacing)));
end

function Rv = rv(prof, coilspacing)
    Rv = 1./((4.*(prof./coilspacing).^2+1.^(1/2)));
end





















