clear, close, clc, 

tic 

% //////////////////////////////////////////
% FORWARD MODEL
% //////////////////////////////////////////

% initial parameters
sig = [1e-3; 5e-2; 3e-2];
nlayer = 100;
depth = 15;
prof = linspace(0, depth, nlayer);

z = [0 ceil([1:length(sig)-1]/length(sig)*depth)];

resmap = ones(length(prof),1)*sig(1);
for i = 1:length(sig)-1
    resmap(round(i/length(sig)*length(resmap)):...
        round((i+1)/length(sig)*length(resmap))) = sig(i+1);
end

%  /////////////////////////////////////////
% sigma_a and data_obs
% ///////////////////////////////////////

% forwardEM1D(ztop, sigma, sep)
% ztop is a column vector containing the z-coordinate of the *top* of each layer [m]
% sigma is a column vector containing the true electrical conductivity of each layer [S/m]
% orient is a scalar value specifying the dipole orientation (0 = vertical; 1 = horizontal)
% sep is a scalar value specifying the coil separation [m]

% valeurs test :
sig = [20e-3; 2e-3;20e-3];
z = [0; 0.5; 1.5];

coilsep = 1:0.5:8; % setting up coilspacing
% orientation = ; % setting up orientation
DataObs = forwardEM1D(z, sig, coilsep); % calling function
save DataObs

% /////////////////////////////////////////////
toc

% /////////////////////////////////
% functions 
% /////////////////////////////////

function forward = forwardEM1D(ztop, sigma, coilsep)
    sep = repmat(coilsep,1, 2)'; % setting up coilspacing
    disp(sep)
    orientation = [ones(1, length(coilsep)) zeros(1, length(coilsep))]; % giving for each coilspacing orientation 1 and 0
    forward = ones(length(sep), 3);
    for i = 1:length(sep)
        if orientation(i) == 1     % 1 = horizontal
            w = ((4.*((ztop./sep(i)).^2)+1).^(1/2))-2.*(ztop./sep(i));
        elseif orientation(i) == 0 % 0 = vertical
            w = ((4.*((ztop./sep(i)).^2)+1).^(1/2)).\1;
        end
        R = w(1:end-1)-w(2:end); % R = weight (R(zn) - R(zn+1))
        R(end+1) = w(end); % R(end) = only R(z)
        surf_sig = sum(sigma(:).*R(:)); % sigma_a
        forward(i,:) = [surf_sig sep(i) orientation(i)]; % building matrix Nx3
    end
end