clear, close, clc, 

% initial parameters
true_sigma = [1/20e3; 1/500e3; 1/20e3]; % 1/[kOhm] -> Hauck et al
true_z = [0; 1; 3];
coilsep = [0.5; 0.5; 1; 1; 2; 2; 4; 4]; % setting up coilspacing
ori =  repmat([0; 1], length(coilsep)/2, 1) ;
data = forwardEM1D(true_sigma, true_z, ori, coilsep); % generate datas
% plot(true_sigma, true_z)
% save datas
save TrueModel true_sigma true_z
save data1D data
