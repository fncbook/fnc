function [mtx, X, Y, vec, unvec, is_boundary] = tensorgrid(x, y)
% TENSORGRID Tensor-product grid.
% Input:
%   x, y           1d projections of the grid nodes (lengths m and n)
% Output:
%   mtx            evaluate a function on the grid (function)
%   X, Y           mtx applied to the coordinate functions (m×n)
%   vec            convert grid shape to vector shape (function)
%   unvec          convert vector shape to grid shape (function)
%   is_boundary    indicator of boundary nodes (logical m×n)
    m = length(x) - 1;
    n = length(y) - 1;
    vec = @(U) U(:);
    unvec = @(u) reshape(u, m+1, n+1);
    [X, Y] = ndgrid(x, y);
    function F = grideval(f)
        F = zeros(size(X));
        for k = 1:numel(X)
            F(k) = f(X(k), Y(k));
        end
    end
    mtx = @grideval;
    
    % Identify boundary points.
    is_boundary = true(m+1, n+1);
    is_boundary(2:m, 2:n) = false;
end