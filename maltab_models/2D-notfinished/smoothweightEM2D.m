function [Dx,Dz] = smoothweightEM2D(nx,nz)
%
% This function computes matrices Wx and Wz that approximate the horizontal
% and vertical second derivatives when applied to a model vector m, where m
% represents a 2D distribution of some geophysical parameter.  Entries of m
% are taken row by row from the 2D image.
%
% Syntax:  [Dx,Dz] = smoothweight(nx,nz);
%
% nx = number of model cells in x direction
% nz = number of model cells in z direction
%
% by James Irving
% June 2019

ncells = nx*nz; % nuber of cell
weights = repmat([1 -2 1],ncells,1); % preset the matrix

% vertical smoothness
Dz = spdiags(weights,[-1 0 1],ncells,ncells);
Dz(1:nz:end,:) = 0; % = 0 avoiding artefacts
Dz(nz:nz:end,:) = 0;

% horizontal smoothness
Dx = spdiags(weights,[-nz 0 nz],ncells,ncells);
Dx(1:nz,:) = 0; % = 0 avoiding artefacts
Dx(((nx-1)*nz+1):nx*nz,:) = 0;