% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

% L is the number of model parameter
% N is the number of datas

tic
clc

disp('Code started')

% load the dataset
load data1D
load TrueModel

% model parameters
ndata = size(data,1); % number of data
nz = 0:.5:20; % vertical discretization
d = data(:, 1); % rename data

nlay = length(nz); % number of discretization layer

% add Gaussian noise (JI: I changed this a bit to have some options)
nperc = 2.5;  % noise level in percent

rng(99999); % set random number seed to have consistent noise
nstd = (nperc/100)*abs(d);                          % noise with a variable standard deviation equal to X% of each data value
% nstd = (nperc/100)*repmat(mean(abs(d)),ndata,1);   % noise with a constant standard deviation equal to X% of mean data value
noise = nstd.*randn(size(d));
d = d + noise;
data(:,1) = d;

% set up m0 and data weight matrix
m0 = true_sigma(1)*ones(size(nz))'; % set up reference model
L = spdiags(1./nstd,0,ndata,ndata); %JI: added data weighting matrix to account for noise characteristics (see inversion notes)
Wd = L'*L;  % setting up Wd

% smoothness
sm = ones(nlay,1);
D = spdiags([sm -2*sm sm], -1:1, nlay, nlay); % build second derivative matrix
D(1, :) = 0; % D(end, :) = 0;

% regularization
alphaz = 1; % setting up the weight of smmoothness and smallness parameters
alphas = .1;
Wm = alphas*speye(length(m0))+ alphaz*(D'*D);

% choose lambda range
lamb = logspace(-2, 12, 2e2); % choose lambda
% lamb = 1e6 ;

% allocate memory for variables
chi2_tot = zeros(size(lamb));
R1D_tot = zeros(size(lamb));
m_tot = cell(size(lamb));

% % inversion iteration
for s = 1:length(lamb)
    [m, G, A] = inversionEM1D(nz, data, lamb(s), Wm, Wd, m0, 1e-10, 2e2); % run inversion
    chi2 = sum(((G*m-d)./nstd).^2)./length(d); % chi-square calculus
    R1d = norm(D*m).^2; % 1D roughness calculus
    
    % store datas
    chi2_tot(:,s) = chi2;
    R1D_tot(:,s) = R1d;
    m_tot{s} = m;
end

% evaluating model uncertainties
disp('Calculating pseudoinverse of A matrix...');

Ainv = pinv(full(A), 1e-10); % invers A
disp('A matrix calculated.');
Cm = Ainv;                            % posterior covariance matrix
R = diag(Ainv*G'*Wd*G);               % model resolution matrix

% adding an arbitrary constrain
ac = 1.67;

% Lagrange parameter selection (chi-square ~= 1, in our case 1.67)
dchi2 = abs(ac - chi2_tot);
ilambda = find(dchi2==min(dchi2));
lambda = lamb(ilambda);

% disp used lambda and chi-square
disp(['optimized lambda = ' num2str(lambda)])
disp(['optimized chi2 = ' num2str(chi2_tot(ilambda))])

% call the model with the best lambda
m = m_tot{ilambda};

% test the inversion by using the forward function
% with the same parameters
data_inv = forwardEM1D(m, nz', data(:,3), data(:,2));

% call and shape true model
true_sigma = [true_sigma(1); true_sigma];
true_z = [true_z; 20];

% bonjour

h_data = data(data(:,3)==1); % select horizontal dipole data
v_data = data(data(:,3)==0); % select vertical dipole data

h_inv_data = data_inv(data_inv(:, 3)==1); % select horizontal predicted data
v_inv_data = data_inv(data_inv(:, 3) == 0); % select vertical predicted data

fig = 1; % remove it if necessary
fig = 1+fig;
invers = figure(fig); clf;
invers.Position = [100 100 1000 600]; % figure size

% m(m(:)<0) = 1e-7;

%----------------------------------------------------
% plot conductivity model dans the true model
invmod = subplot(2,2,1);
stairs(true_sigma, true_z, '--k')
hold on
stairs(m,nz, 'b')
title('(a)')% Subsurface conductivity model')
xlabel('\sigma [Sm^{-1}]')
ylabel('Depth [m]')
ylim([-0.3 10])
xlim padded
set(gca,'Ydir','reverse')
set(gca,'XScale','log')
grid on
legend('true model', 'predicted', 'Location', 'southwest')

%---------------------------------------------------------
% plot L-curve
subplot(2,2,4)
loglog(chi2_tot, R1D_tot, 'b.', 'MarkerSize', 1)
hold on
loglog(chi2_tot(ilambda), R1D_tot(ilambda), 'o', 'MarkerSize', 11)
hold off
title('(e)')
xlabel('\chi^2')
ylabel('Roughness')
legend('L - curve', ['lambda ', num2str(lambda, '%.e')], 'Location', 'best')
axis padded
grid on

%---------------------------------------------------------
% plot model resolution
Rmod = subplot(2, 2, 2);
semilogx(R, nz, 'b')
title('(d)')
xlabel('Resolution')
ylabel('Depth [m]')
% axis padded
xlim([0 1])
ylim([-0.3 10])
set(gca,'Ydir','reverse')
grid on 

%----------------------------------
% plot vertical true data and predicted data
nd = unique(data(:, 2));
dav = subplot(2,4,5);
plot(nd, v_data, 'kv');
hold on
plot(nd, v_inv_data, 'bx')
title('(b)')
legend('observed', 'predicted', 'Location', 'northeast')
xlabel('Coilspacing [m]')
ylabel('\sigma_a [Sm^{-1}]')
% xlim([0 0.025])
set(gca,'YScale','log')
axis padded
grid on

% plot horizontal data and predicted data
dah = subplot(2,4,6);
plot(nd, h_data, 'ko');
hold on
plot(nd, h_inv_data, 'bx')
title('(c)')
legend('observed', 'predicted', 'Location', 'northeast')
xlabel('Coilspacing [m]')
set(gca,'YScale','log')
axis padded
grid on

linkaxes([dav dah], 'xy')

%-------------------------------------

linkaxes([invmod, Rmod], 'y')

% give a big title
sgtitle(['\lambda = ', num2str(lamb(ilambda), '%.e'), ',  \chi^2 = ',...
    num2str(chi2_tot(ilambda)), ',  p_G(d) = ', num2str(nperc), ' %', ',  \alpha_s = ', num2str(alphas), ',  \alpha_z = ', num2str(alphaz)], 'FontSize', 12)

% save figure
% saveas(invers, 'inv-arte-1d.eps', 'epsc')
% saveas(invers, 'inv-arte-1d.fig')

disp('code finished : ')
toc