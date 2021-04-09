clear, close, clc, 

tic
% initial parameters
sig = [20e-3; 2e-3; 10e-3];
z = [0; 0.5; 1.5];
coilsep = 0.1:0.1:8; % setting up coilspacing
data = forwardEM1D(sig, z, coilsep); % generate datas

% save datas
save InitialModel sig
save data1D data
toc