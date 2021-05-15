% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

% L is the number of model parameter
% N is the number of datas

tic
clear; close all; clc; tic

% load the dataset
load data1D
load TrueModel

% model parameters
ndata = size(data,1);
nz = 0:.2:10;
nlay = length(nz);
d = data(:, 1);

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
alphas = 1;
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

% inversion normal
% [m, G, A] = inversionEM1D(nz, data, lamb, Wm, Wd, m0, 1e-10, 2e2);
% (chi2<1, fitting data too well; chi2>1, not fitting data enough)
% chi2 = sum(((G*m-d)./nstd).^2)./length(d);
% disp(['Chi-squared misfit statistic = ',num2str(chi2)]);
% disp(['lambda = ',num2str(lamb, '%.e')])

% evaluating model uncertainties
disp('Calculating pseudoinverse of A matrix...');
Ainv = pinv(full(A), 1e-10);
disp('A matrix calculated.');
Cm = diag(Ainv);                      % posterior covariance matrix
R = diag(Ainv*G'*Wd*G);               % model resolution matrix
% R_mat = G'*(G*G')^-1 *G;               % model resolution matrix
% R = diag(R);

% Lagrange parameter selection
dchi2 = abs(chi2_tot-ndata);
ilambda = find(dchi2==min(dchi2));
lambda = lamb(ilambda);

% disp used lambda and chi-square
disp(['optimized lambda = ' num2str(lambda)])
disp(['optimized chi2 = ' num2str(chi2_tot(ilambda))])

% call the model with the best lambda
m = m_tot{ilambda};

% call and shape true model
true_sigma = [true_sigma(1); true_sigma];
true_z = [true_z; 20];

figure(1)
subplot(2,2,2)
scatter(chi2_tot-ndata, R1D_tot,12, '.')
hold on
scatter(chi2_tot(ilambda)-ndata, R1D_tot(ilambda), 100, 'o')
hold off
title('L-curve')
xlabel('\chi^2 - N')
ylabel('R1D')
legend('L - curve', ['lambda ', num2str(lamb(ilambda), '%.e')])
axis padded
grid on

% plot
% invers = figure(2);
subplot(2,2,1)
stairs(m,nz)
hold on
stairs(true_sigma, true_z)
title('Subsurface conductivity model')
subtitle(['\lambda = ', num2str(lamb(ilambda), '%.e'), ',  \chi^2 = ',...
    num2str(chi2_tot(ilambda)), ',  p_G(d) = ', num2str(nperc), ' %'])
xlabel('conductivity \sigma_a [Sm^{-1}]')
ylabel('depth [m]')
legend('modeled', 'observed', 'location', 'SouthWest')
ylim([-1 6])
xlim padded
% set(gca,'XScale','log')
xticks([0 0.005 0.01 0.015 0.02 0.025 0.03])
set(gca,'Ydir','reverse')
grid on
% saveas(invers, 'inversion.png')

subplot(2, 2, 3)
semilogy(R)
title('model resolution')
ylabel('resolution rate')
xlabel('model parameters L')
axis padded
grid on

subplot(2, 2, 4)
semilogy(Cm)
title('covariance model')
ylabel('covariance rate')
xlabel('model parameters L')
axis padded
grid on

disp('code finished : ')
toc