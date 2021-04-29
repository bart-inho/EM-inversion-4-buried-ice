clear, close, clc, 

tic
% initial parameters
xlog = 0:0.5:20; %[m] horizontal discretization
nmeasure = length(xlog); % number of horizontal measurments
ztop = repmat([0; 0.5; 3; 7], 1, nmeasure); % top layer vertical coordinate

sig = repmat([20e-3; 1e-3; 20e-3; 10e-3], 1, nmeasure); % true model map

coilsep = repmat(0.1:0.1:10, nmeasure, 1)'; % coilseparations
ori = repmat([0 1], length(xlog), size(coilsep, 1)/2)'; % orientation of the dipole (0 = vertical, 1 = horizontal)

data = []; % preset data matrix

for i = 1:length(xlog)
    % generate datas in a matrix that contains physical properties
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))];
end

% save datas
save data2D data
toc