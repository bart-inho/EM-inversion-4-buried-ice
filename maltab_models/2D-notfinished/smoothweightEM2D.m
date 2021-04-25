function [Dz,Dx] = smoothweightEM2D(nx,nz)
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

ncells = nx*nz;
weights = repmat([1 -2 1],ncells,1);
Dz = spdiags(weights,[-1 0 1],ncells,ncells);
Dz(1:nx:end,:) = 0;
Dz(nx:nx:end,:) = 0;
Dx = spdiags(weights,[-nx 0 nx],ncells,ncells);
Dx(1:nx,:) = 0;
Dx(((nz-1)*nx+1):nx*nz,:) = 0;