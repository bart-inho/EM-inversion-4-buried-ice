% resolution, phi and covariance matrix matrix
tic
% load the dataset

load data2D
noi = 0.025; % *100 = n%

xlog = unique(data(:,4)); % horizontal size of the model
ztop = 0:0.05:5; % vertical size of the model
nx = length(xlog); % number of discretization layers horizontal
nz = length(ztop); % number of discretization layers vertical



% M = 1./(G'*G+lamb*Wm)*G';
% % I = diag(ones(nx*nz, 1));
% 
% cov_mat = noi^2 *M*M';

cov_mat = cov(m);

figure()
imagesc(cov_mat)
toc