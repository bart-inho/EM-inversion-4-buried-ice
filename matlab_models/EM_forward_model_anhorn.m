clear all, close all, clc, 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////
tic

% initial conductivity
sig = [1e-3; 1e-4; 6e-4];
depth = 10;

% frequency and depth
freq = 1.6e5;
% depth = round(500/sqrt(sig(1)*freq)); % depth [m]
z = [];
for i = 1:length(sig)
    z = [z ceil(i/length(sig) * depth)];
end

% building inital conductivity model
distance = 50; % A-B distance [m]
resmap = ones(depth, distance);
% resmap(1:end,1:end) = sig(1);
% resmap(z(1):z(2),1:end) = sig(2);
% resmap(z(2):end,1:end) = sig(3);
resmap(1:end,1:end) = sig(1);
for i = 1:length(z)-1
    resmap(z(i):z(i+1), 1:end) = sig(i+1);
end

% plotting the initial model 
% initial_model = figure(1);
% hold on
% pcolor(resmap)
% xlabel('x')
% ylabel('depth')
% title('conductvity model')
% colorbar
% hold off
% saveas(initial_model, 'initial_model.png')

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

%  /////////////////////////////////////////
% giving sigma_a
% ///////////////////////////////////////

% forwardEM1D(ztop, sigma, orient, sep)
% ztop is a column vector containing the z-coordinate of the *top* of each layer [m]
% sigma is a column vector containing the true electrical conductivity of each layer [S/m]
% orient is a scalar value specifying the dipole orientation (0 = vertical; 1 = horizontal)
% sep is a scalar value specifying the coil separation [m]
% if two direction desired : sigma_a = [sigma_a_v sigma_a_h] ?

coilsep = 20e1;
ori = 0;
sigma_a = forwardEM1D(z, sig, ori, coilsep);
disp(sigma_a)

% /////////////////////////////////////////////
% Inversion ???
% /////////////////////////////////////////////

dstd = 1e0; % 
lambda = 1e0; % 
k = 1:1:length(centroid);
m0 = ones(length(centroid), 1) * mean(sig);

% not sure of how to build matrix A (or G depending the source)
A = zeros(length(centroid), length(k));
A(:, :) = centroid(:).^repmat(k, length(centroid), 1);

L = (1/dstd)*speye(length(centroid));
Qv = (4*(centroid./coilsep))./((4*(centroid.^2)./(coilsep.^2)+1).^(3/2));
Qh = 2 - ((4*(centroid./coilsep))/(sqrt(4*((centroid.^2)./(coilsep.^2))))) ;

Wd = L'*L;
Wm = diag(1./Qv);

% We get a correct matrix m^N but not optimized.
m = m0 + inv(Wm)*A'*inv(A*inv(Wm)*A'+Wd*lambda^2)*(sigma_a-A*m0);
% m(m<0) = 0; % to not have "negative conductivity ?

% optimization ?
% phid = (A*m - sigma_a)'.*Wd*(A*m - sigma_a);
% phim = (m - m0)'.*Wm*(m-m0);

inversion = figure(4);
hold on
plot(m, centroid, 'r')
plot(m0, centroid, 'k-.')
plot(resEM, centroid, 'b--')
legend('inverted model', 'm0', 'initial model')
ylabel('centroid')
xlabel('\sigma [S/m]')
% xlim([-1e-4 2e-3])
hold off
saveas(inversion, 'inversion_plt.png')

% ??

% /////////////////////////////////////////////
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