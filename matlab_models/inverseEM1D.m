% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear, close all, clc, tic

% load the dataset
load data1D
load InitialModel

% model parameters
ndata = size(data,1);
ztop = 0:0.1:2;
nlay = length(ztop);
sigma_a = data(:, 1);

figure()
plot(sigma_a, '*')
hold on

% noise
for i = 1:length(sigma_a)
    sigma_a(i) = sigma_a(i)+normrnd(0,1)*5e-5;
end

plot(sigma_a, 'x')
hold off

% G matrix
G = zeros(ndata, nlay);   % initialize W matrix
for i=1:size(G, 1)        % populate the matrix row-by-row
    R = Rcalc(ztop, data(i,2), data(i,3));
    G(i, :) = R;
end

% smallness
I = eye(nlay);

% smoothness
n = nlay;           % the size of the matrix 
v = ones(n,1)*(-2);   % make the vector for the main diagonal
u = ones(n-1,1); % make the vector for +1 and -1 diagonal
D = diag(v)+diag(u,1)+diag(u,-1); % combine everything together

% regularization parameters
% alpha_s = 1 ;
% alpha_z = 0.1 ;

% regularization
% Wm = alpha_s*I + alpha_z*(D'*D);
% lamb = 1e-3;

% mark
mark = 1;

% choose the model
if mark == 1
    Wm = D'*D;
    lamb = 1e-4;
elseif mark == 0
    Wm = I;
    lamb = 1e-3;
end

% setting up inverse model
A = G'*G + lamb*Wm;
b = G'*sigma_a;

% inversion
m = cgs(A, b, 1e-10, 5e2);

% initial model we know
sigma = ones(length(ztop), 1)*sig(1);
sigma(ceil(1/3*length(ztop)):ceil(2/3*length(ztop))) = sig(2);
sigma(ceil(2/3*length(ztop)):end) = sig(3);
toc

% plot
invers = figure();
plot(m,ztop)
hold on
plot(sigma, ztop)
legend('inverted model', 'initial model', 'location', 'best')
set(gca,'Ydir','reverse')
saveas(invers, 'inversion.png')