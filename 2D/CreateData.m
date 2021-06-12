clear; close all; clc;

tic
% initial parameters
disc = 0.5;
xlog = 0:disc:40; %[m] horizontal discretization
nmeasure = length(xlog); % number of horizontal measurments

% bloc depth
z1bloc1 = 4 + disc;
z2bloc1 = 9;

z1bloc2 = 1+ disc;
z2bloc2 = 5;

% 1 bloc
% ztop = repmat([0; 1; 3], 1, nmeasure); % top layer vertical coordinate

% two blocs
ztop = [repmat([0; z1bloc1/disc; z2bloc1/disc], 1, ceil(nmeasure/2)) repmat([0; z1bloc2/disc; z2bloc2/disc], 1, round(nmeasure/2))]; % top layer vertical coordinate
mo_sig = 1/20e3;
ic_sig = 1/500e3;

% bloc shape coordinate
x1bloc1 = 18;
x2bloc1 = 24;
x1bloc2 = 52;
x2bloc2 = 62;

% 2 blocs
sig = mo_sig*ones(3, nmeasure);
sig(2:2, x1bloc1:x2bloc1) = ic_sig;
sig(2:2, x1bloc2:x2bloc2) = ic_sig;

% 1 bloc
% sig(2:2,39:43) = ic_sig;

coilsep = repmat([0.5 0.5 1 1 2 2 4 4], nmeasure, 1)'; % coilseparations
ori = repmat([1 0], nmeasure, 4)'; % orientation of the dipole (0 = vertical, 1 = horizontal)
imagesc(sig)
data = []; % preset data matrix

for i = 1:nmeasure
    % generate datas in a matrix that contains physical properties
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))];
end

zdis = 0:.5:20;
ndis = length(zdis);

% 2 blocs
sigini = mo_sig*ones(ndis, nmeasure);
sigini(z1bloc1/disc:z2bloc1/disc, x1bloc1:x2bloc1) = ic_sig;
sigini(z1bloc2/disc:z2bloc2/disc, x1bloc2:x2bloc2) = ic_sig;
 
% 1 bloc
% sigini = mo_sig*ones(ndis, nmeasure);
% sigini(2:6,39:43) = ic_sig;

pcolor(xlog, zdis, sigini)
set(gca,'Ydir','reverse')
axis image

% save datas
save zdis zdis
save sigini sigini
save data2D data
toc