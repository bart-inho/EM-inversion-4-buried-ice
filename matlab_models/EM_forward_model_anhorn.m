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
plot(resEM(1:end), centroid(1:end))
xlim([0,1.2e-3])
xlabel('Conductivity \sigma')
ylabel('centroid')
title('forward model')

% model for later -> no indexing inversion is necessary, resEM is already ok.
model = thk.*resEM;

%/////////////////////////////////////////
% giving sigma_a
% ///////////////////////////////////////

coilspace = [5 8 10 12 14]; % coil spacing [m]
z1 = 10; z2 = 20; % depth [m]

sigma_a_h = sig1*(1-rh(z1, coilspace)) + sig2*(rh(z1, coilspace) - rh(z2, coilspace)) + sig3*rh(z2, coilspace);
disp('horizontal sigma a = '); disp(sigma_a_h);
sigma_a_v = sig1*(1-rv(z1, coilspace)) + sig2*(rv(z1, coilspace) - rv(z2, coilspace)) + sig3*rv(z2, coilspace);
disp('vertical sigma a = '); disp(sigma_a_v);

% functions 
function Rh = rh(prof, coilspacing)
    Rh = 1./((4.*(prof./coilspacing).^2+1.^(1/2)-2.*(prof./coilspacing)));
end

function Rv = rv(prof, coilspacing)
    Rv = 1./((4.*(prof./coilspacing).^2+1.^(1/2)));
end