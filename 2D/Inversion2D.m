% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clc; tic, close all
% load the dataset
load data2D
load sigini
load zdis

xlog = unique(data(:,4)); % horizontal size of the model
% zdis = 0:.5:20; % vertical size of the model
nx = length(xlog); % number of discretization layers horizontal
nz = length(zdis); % number of discretization layers vertical
ndata = size(data,1); % number of sigma_a

% add Gaussian noise
nperc = 2.5;  % noise level in percent
rng(99999); % set random number seed to have consistent noise
d = data(:,1); % apparent conductivity data
nstd = (nperc/100)*abs(d);                          % noise with a variable standard deviation equal to X% of each data value
%nstd = (nperc/100)*repmat(mean(abs(d)),ndata,1);   % noise with a constant standard deviation equal to X% of mean data value
noise = nstd.*randn(size(d));
d = d + noise;
data(:,1) = d;

% set parameters for the regularization
% alphax should generally be much bigger than alphaz for layered media
alphas = .1;  % weight on model smallness relative to reference model (see inversion notes)
m0 = 1/20e3*ones(nx*nz,1); % reference constant conductivity model (not considered if alphas=0)
alphax = 1; % weight on model smoothness in x-direction
alphaz = 10;  % weight on model smoothness in z-direction
% set reference model

% calculate data and model weighting matrices
L = spdiags(1./nstd,0,ndata,ndata); %JI: added data weighting matrix to account for noise characteristics (see inversion notes)
Wd = L'*L;  
[Dx,Dz] = smoothweightEM2D(nx,nz);
Wm = alphas*speye(length(m0)) + alphax*(Dx'*Dx) + alphaz*(Dz'*Dz);

% inversion parameters knowing that lambda
lamb = logspace(3, 15, 1e2); % trade-off parameter
% lamb = 10;
tol = 1e-9; % tolerance for conjugate gradient solver
itmax = 500; % maximum # of iterations

% allocate memory to variables
chi2_tot = zeros(size(lamb));
R1D_tot = zeros(size(lamb));
m_tot = cell(size(lamb));
A_tot = cell(size(lamb));
G_tot = cell(size(lamb));

for s = 1:length(lamb)
    [m, m1, G, A] = inversionEM2D(zdis, nz, nx, data, lamb(s), Wm, Wd, m0, tol, itmax);
    
    % chi2
    chi2 = sum(((G*m1-d)./nstd).^2)/length(d); % chi-square
    R = norm(Wm*m1).^2; % resolution matrix
    
    % store datas
    chi2_tot(:,s) = chi2;
    R1D_tot(:, s) = diag(R);
    m_tot{s} = m;
    A_tot{s} = A;
    G_tot{s} = G;
end

% Find best lambda index
b_chi2value = 1;
dchi2 = abs(b_chi2value - chi2_tot);
ilambda = find(dchi2==min(dchi2));
lambda = lamb(ilambda);

% display important values
disp(['Choosen lambda = ', num2str(lambda,'%.e')])
disp(['best chi2 value = ' num2str(b_chi2value)])
disp(['Chi-squared misfit statistic = ',num2str(chi2_tot(:, ilambda))]);

% select the best model
m = m_tot{ilambda};
A = A_tot{ilambda};
G = G_tot{ilambda};

% evaluating model uncertainties
disp('Calculating pseudoinverse of A matrix..., please wait');
Ainv = pinv(full(A),1e-10);
disp('A matrix calculated.');
% Cm = Ainv;                      % posterior covariance matrix
% dCm = reshape(diag(Cm),nz,nx);
R = Ainv*G'*Wd*G;               % model resolution matrix
dR = reshape(diag(R),nz,nx);

% plot
invers = figure(1);
invers.Position = [100 100 1000 600];
subplot(2, 1, 1)
tiledlayout(2,1);
nexttile
inv = pcolor(xlog, zdis(1:21), m(1:21,:));
inv.EdgeColor = 'none';
title('(a)')
xlabel('[m]')
ylabel('depth [m]')
axis image
caxis([1/500e3 1/20e3])
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = '\sigma [S/m]';
set(gca, 'YDir','reverse')
nexttile
true = pcolor(xlog, zdis(1:21),sigini(1:21,:));
set(gca, 'YDir','reverse')
true.EdgeColor = 'none';
title('(b)')
xlabel('[m]')
ylabel('depth [m]')
axis image

sgtitle(['\lambda = ', num2str(lambda, '%.e'),...
    ',  \chi^2 = ', num2str(chi2_tot(ilambda)), ',  P_G(m) = ', num2str(nperc), ' %'],...
    'FontSize', 12)

resolu = figure(2);
resolu.Position = [100 100 1000 600];
pcolor(xlog,zdis(1:21),dR(1:21, :));
% title('Resolution');
xlabel('Position [m]')
ylabel('Depth [m]')
axis image
c = colorbar;
c.Label.String = 'Resolution';
caxis([1e-4 1])
set(gca,'ColorScale','log')
true.EdgeColor = 'none';
set(gca, 'YDir','reverse')

lcurve = figure(3);
lcurve.Position = [100 100 1000 600];
loglog(chi2_tot, R1D_tot,'b.', 'Markersize', 5)
hold on
loglog(chi2_tot(ilambda), R1D_tot(ilambda), 'ro', 'Markersize', 10)
axis equal
% title('L-curve')
legend('L-curve', ['lambda = ' num2str(lambda, '%.e')])
xlabel('\chi^2')
ylabel('roughness R')

% save figures
% saveas(invers, 'inv-seminoise-2d.eps', 'epsc')
% saveas(invers, 'inv-seminoise-2d.fig')

disp('code finished :')
toc