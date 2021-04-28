clear, close, clc, 

tic
% initial parameters
xlog = 0:0.5:25; %[m]
nlay = length(xlog);
sig = repmat([20e-3; 1e-3; 20e-3; 10e-3], 1, nlay);
ztop = repmat([0; 1; 2; 8], 1, nlay);
coilsep = repmat(0.1:0.1:10, nlay, 1)';
ori = repmat([0 1], length(xlog), size(coilsep, 1)/2)';

data = [];

for i = 1:length(xlog)
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))]; % generate datas
end

% save datas
save data2D data
toc