clear, close, clc, 

tic
% initial parameters
xlog = 0:1:20; %[m] horizontal discretization
nmeasure = length(xlog); % number of horizontal measurments
ztop = repmat([0; 1.5; 3], 1, nmeasure); % top layer vertical coordinate

sig = repmat([0.02; 0.001; 0.02], 1, nmeasure); % true model map
coilsep = repmat([0.5 0.5 1 1 2 2 4 4], nmeasure, 1)'; % coilseparations
ori = repmat([1 0], nmeasure, 4)'; % orientation of the dipole (0 = vertical, 1 = horizontal)
imagesc(sig)
data = []; % preset data matrix

for i = 1:nmeasure
    % generate datas in a matrix that contains physical properties
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))];
end

% save datas
save data2D data
toc