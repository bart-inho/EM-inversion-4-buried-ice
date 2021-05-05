clear, close, clc, 

tic
% initial parameters
xlog = 0:0.2:10; %[m] horizontal discretization
nmeasure = length(xlog); % number of horizontal measurments
ztop = repmat([0; 1; 2], 1, nmeasure); % top layer vertical coordinate

sig = repmat([0.02; 0.002; 0.02], 1, nmeasure); % true model map
coilsep = repmat([4 4 6 6 8 8], nmeasure, 1)'; % coilseparations
ori = repmat([1 0], nmeasure, 3)'; % orientation of the dipole (0 = vertical, 1 = horizontal)

data = []; % preset data matrix

for i = 1:nmeasure
    % generate datas in a matrix that contains physical properties
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))];
end

% save datas
save data2D data
toc