% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; clc; close all
% load the dataset
load data2D

xlog = unique(data(:,4)); % horizontal size of the model
ztop = 0:0.05:5; % vertical size of the model
nx = length(xlog); % number of discretization layers horizontal
nz = length(ztop); % number of discretization layers vertical
ndata = size(data,1); % number of sigma_a

% add Gaussian noise (JI: I changed this a bit to have some options)
nperc = 1;  % noise level in percent
rng(99999); % set random number seed to have consistent noise
d = data(:,1); % apparent conductivity data
nstd = (nperc/100)*abs(d);                          % noise with a variable standard deviation equal to X% of each data value
%nstd = (nperc/100)*repmat(mean(abs(d)),ndata,1);    % noise with a constant standard deviation equal to X% of mean data value
noise = nstd.*randn(size(d));
d = d + noise;

% setting up G matrix
ndata_p = size(data(data(:, 4) == 1), 1); % select data for each measurment point
G_stack = cell(1, nz); % predefine the list that will contain all 1D G marix

for j = 1:nx
    G_prov = zeros(ndata_p, nz);   % initialize W matrix
    for i=1:size(G_prov, 1)        % populate the matrix row-by-row
        R = weightEM2D(ztop, data(i,2), data(i,3));
        G_prov(i, :) = R;
    end
    G_stack{j} = G_prov;
end

% G matrix
G = sparse(blkdiag(G_stack{:})); % stack, blockdiag and sparse

% set parameters for the regularization
% JI: alphax should generally be much bigger than alphaz for layered media
alphas = 0;  % weight on model smallness relative to reference model (see inversion notes)
m0 = 0.02*ones(nx*nz,1); % reference constant conductivity model (not considered if alphas=0)
alphax = 10; % weight on model smoothness in x-direction
alphaz = 1; % weight on model smoothness in z-direction
% set reference model

% calculate data and model weighting matrices
L = spdiags(1./nstd,0,ndata,ndata); %JI: added data weighting matrix to account for noise characteristics (see inversion notes)
Wd = L'*L;  
[Dx,Dz] = smoothweightEM2D(nx,nz);
Wm = alphas*speye(length(m0)) + alphax*(Dx'*Dx) + alphaz*(Dz'*Dz);

% inversion parameters
lamb = 1e6; % trade-off parameter
tol = 1e-10; % tolerance for conjugate gradient solver
itmax = 500; % maximum # of iterations

% solve linear system  %JI: note additions here to reflect above changes
A = G'*Wd*G + lamb*Wm; 
b = G'*Wd*d + lamb*Wm*m0;
[m1,fl1,rr1,it1,rv1] = cgs(A, b, tol, itmax);
m = reshape(m1,nz,nx);

% calculate chi-squared misfit statistic (JI: added this part...
% (chi2<1, fitting data too well; chi2>1, not fitting data enough)
chi2 = sum(((G*m1-d)./nstd).^2)./length(d);
disp(['Chi-squared misfit statistic = ',num2str(chi2)]);

% evaluating model uncertainties
disp('Calculating pseudoinverse of A matrix...');
Ainv = pinv(full(A),1e-10);
Cm = Ainv;                      % posterior covariance matrix
dCm = reshape(diag(Cm),nz,nx);
R = Ainv*G'*Wd*G;               % model resolution matrix
dR = reshape(diag(R),nz,nx);


% figure()
% semilogy(0:length(rv1)-1,rv1/norm(b),'-')
% hold on
% yline(tol,'r--');
% xlabel('Iteration number')
% ylabel('Relative residual')

figure
subplot(3,1,1)
%inv = pcolor(xlog, ztop, m);
imagesc(xlog,ztop,m);
title('Inverted model')
%subtitle({['\lambda = ' num2str(lamb) ', \alpha_x = ' num2str(alphax), ', \alpha_z = ' num2str(alphaz)]...
%    ['nx = ' num2str(nx) ', nz = ' num2str(nz) ', %noise = ' num2str(nperc)]...
%    ['coilspacing = [ ' num2str(data(1, 2)) ', ' num2str(data(3, 2)) ', ' num2str(data(5, 2)) ' ]']})
set(gca, 'YDir','reverse')
xlabel('Position [m]')
ylabel('Depth [m]')
%shading interp;
%axis image
c = colorbar;
c.Label.String = '\sigma [S/m]';
% caxis([0 max(m1)]);
%saveas(inv, 'inversion_2d', 'png')

subplot(3,1,2)
imagesc(xlog,ztop,log10(dCm));
title('log10(Posterior parameter variance)');
xlabel('Position [m]')
ylabel('Depth [m]')
colorbar

subplot(3,1,3)
imagesc(xlog,ztop,log10(dR));
title('log10(Resolution)');
xlabel('Position [m]')
ylabel('Depth [m]')
colorbar
