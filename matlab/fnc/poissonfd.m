function [X, Y, U] = poissonfd(f,g,m,xspan,n,yspan)
%POISSONFD   Solve Poisson's equation by finite differences.
% Input:
%   f            forcing function (function of x,y)
%   g            boundary condition (function of x,y)
%   m            number of grid points in x (integer)
%   xspan        endpoints of the domain of x (2-vector)
%   n            number of grid points in y (integer)
%   yspan        endpoints of the domain of y (2-vector)
%
% Output:
%   U            solution (m+1 by n+1)
%   X,Y          grid matrices (m+1 by n+1)

% Discretize the domain.
[x, Dx, Dxx] = diffmat2(m, xspan);
[y, Dy, Dyy] = diffmat2(n, yspan);
[mtx, X, Y, vec, unvec, is_boundary] = tensorgrid(x, y);

% Form the collocated PDE as a linear system. 
Ix = speye(m+1);  Iy = speye(n+1);
A = kron(Iy, sparse(Dxx)) + kron(sparse(Dyy), Ix);  % Laplacian matrix
b = vec(mtx(f));

% Replace collocation equations on the boundary.
scale = max(abs(A(n+2, :)));
I = speye(size(A));
idx = vec(is_boundary);
A(idx, :) = scale * I(idx, :);           % Dirichet assignment
b(idx) = scale * g( X(idx),Y(idx) );     % assigned values

% Solve the linear sytem and reshape the output.
u = A \ b;
U = unvec(u);