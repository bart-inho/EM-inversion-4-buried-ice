clear, close, clc, 

tic
% initial parameters
xlog = 0:0.2:10; %[m]
nlay = length(xlog);
sig = repmat([20e-3; 2e-3; 10e-3], 1, nlay);
ztop = repmat([0; 1; 1.5], 1, nlay);
coilsep = repmat(0.1:0.1:8, nlay, 1)';
ori = repmat([0 1], length(xlog), length(coilsep)/2)';

data = [];

for i = 1:length(xlog)
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))]; % generate datas
end

% save datas
save data2D data
toc