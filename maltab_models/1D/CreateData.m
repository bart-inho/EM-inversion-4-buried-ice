clear, close, clc, 

% initial parameters
sig = [20e-3; 2e-3; 10e-3];
z = [0; 0.5; 2];
coilsep = repmat(4:0.5:8, 1, 2)'; % setting up coilspacing
ori =  repmat([0; 0], length(coilsep)/2, 1) ;
data = forwardEM1D(sig, z, ori, coilsep); % generate datas

% save datas
save TrueModel sig
save data1D data