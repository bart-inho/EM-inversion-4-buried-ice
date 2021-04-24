% simple 1D inversion of FDEM data acquired over a layered Earth using
% multiple coil separations and both vertical and horizontal dipoles

clear; close all; clc; tic

% load the dataset
load data2D

xlog = 0:1:10; %

% model parameters
ztop = 0:0.1:2;
sigma_a = data(:, 1);

ndata = size(data, 1);
ndata_pos = size(data(data(:, 4) == 0),1);
nlay = length(ztop);
G = zeros(ndata, nlay*ndata_pos);   % initialize W matrix
G_prov = zeros(ndata_pos, nlay);

for j = 
    for i=1:size(G_prov, 1)        % populate the matrix row-by-row
        R = WeightEM(ztop, data(i,2), data(i,3));
        G_prov(i, :) = R;
    end
    
