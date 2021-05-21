% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

% L is the number of model parameter
% N is the number of datas

tic
clear; close all; clc; tic

disp('Code started')
% load the dataset
load data1D
load TrueModel

% model parameters
ndata = size(data,1);
nz = 0:.25:10;
d = data(:, 1);

% % Creation of the structure in depth
% nlayer = 16; % Number of layers
% 
% % Thicknesses of layers [m]
% thick = ones(nlayer,1);
% thick(1) = 0.1;
% for j=2:nlayer-1
%     thick(j) = 1.2*thick(j-1);
% end, clear j
% thick(end) = 10;
% 
% % Depths of layer interfaces [m]
% nz = zeros(size(thick,1)+1,1);
% for i = 1:length(thick)
%     nz(i+1) = nz(i)+thick(i);
% end
% nz = nz';

nlay = length(nz);

% add Gaussian noise (JI: I changed this a bit to have some options)
nperc = 2.5;  % noise level in percent

rng(99999); % set random number seed to have consistent noise
nstd = (nperc/100)*abs(d);                          % noise with a variable standard deviation equal to X% of each data value
%nstd = (nperc/100)*repmat(mean(abs(d)),ndata,1);   % noise with a constant standard deviation equal to X% of mean data value
noise = nstd.*randn(size(d));
d = d + noise;
data(:,1) = d;

% set up m0 and data weight matrix
m0 = 0.02*ones(size(nz))';
L = spdiags(1./nstd,0,ndata,ndata); %JI: added data weighting matrix to account for noise characteristics (see inversion notes)
Wd = L'*L;  

% smoothness
sm = ones(nlay,1);
D = spdiags([sm -2*sm sm], -1:1, nlay, nlay);
D(1, :) = 0; %D(end, :) = 0;

% regularization
alphaz = 1;
alphas = .1;
Wm = alphas*speye(length(m0))+ alphaz*(D'*D);

% choose lambda range
lamb = logspace(-2, 8, 1e3);
% lamb = 1e6 ;

% allocate datas
chi2_tot = zeros(size(lamb));
R1D_tot = zeros(size(lamb));
m_tot = cell(size(lamb));

% % inversion iteration
for s = 1:length(lamb)
    [m, G, A] = inversionEM1D(nz, data, lamb(s), Wm, Wd, m0, 1e-10, 2e2);
    chi2 = sum(((G*m-d)./nstd).^2)./length(d);
    R = norm(D*m).^2;
    
    % store datas
    chi2_tot(:,s) = chi2;
    R1D_tot(:,s) = R;
    m_tot{s} = m;
end

% evaluating model uncertainties
disp('Calculating pseudoinverse of A matrix...');
Ainv = pinv(full(A), 1e-10);
disp('A matrix calculated.');
Cm = Ainv;                            % posterior covariance matrix
R = diag(Ainv*G'*Wd*G);               % model resolution matrix

% Lagrange parameter selection
dchi2 = abs(chi2_tot-ndata);
ilambda = find(dchi2==min(dchi2));
lambda = lamb(ilambda);

% disp used lambda and chi-square
disp(['optimized lambda = ' num2str(lambda)])
disp(['optimized chi2 = ' num2str(chi2_tot(ilambda))])

% call the model with the best lambda
m = m_tot{ilambda};

% test the inversion by using the forward function
% with the same parameters
sigma_a_inv = forwardEM1D(m, nz', data(:,3), data(:,2));

% call and shape true model
true_sigma = [true_sigma(1); true_sigma];
true_z = [true_z; 20];

invers = figure(1);
invers.Position = [100 100 1000 600];

subplot(2,2,1)
stairs(m,nz)
hold on
stairs(true_sigma, true_z)
title('Subsurface conductivity model')
subtitle(['\lambda = ', num2str(lamb(ilambda), '%.e'), ',  \chi^2 = ',...
    num2str(chi2_tot(ilambda)), ',  p_G(d) = ', num2str(nperc), ' %'])
xlabel('conductivity \sigma [Sm^{-1}]')
ylabel('depth [m]')
legend('modeled', 'observed', 'location', 'SouthWest')
ylim([-0.1 6])
xlim([0 0.025])
% set(gca,'XScale','log')
xticks([0 0.005 0.01 0.015 0.02 0.025 0.03])
set(gca,'Ydir','reverse')
grid on

subplot(2,2,3)
plot(sigma_a_inv(:, 1),1:1:ndata, '-o')
hold on
plot(data(:, 1),1:1:ndata, '-x')
% xticks([0 0.005 0.01 0.015 0.02 0.025 0.03])
title('Apparent conductivity')
legend('synthetic data', 'inverted data')
ylabel('N data')
xlabel('\sigma_a [Sm^{-1}]')
xlim([0 0.025])
ylim padded
grid on

subplot(2,2,2)
scatter(chi2_tot, R1D_tot,1, '.')
hold on
scatter(chi2_tot(ilambda), R1D_tot(ilambda), 100, 'o')
hold off
title('L-curve')
xlabel('\chi^2')
ylabel('Roughness rate')
set(gca,'Xscale','log')
set(gca,'Yscale','log')
legend('L - curve', ['lambda ', num2str(lamb(ilambda), '%.e')])
axis padded
grid on

subplot(2, 2, 4)
semilogy(R)
title('Model resolution')
ylabel('Resolution rate')
xlabel('L model parameters')
axis padded
grid on 

saveas(invers, 'inversion.png')

disp('code finished : ')
toc