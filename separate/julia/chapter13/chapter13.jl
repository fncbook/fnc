# begin tensorgrid
"""
    tensorgrid(x, y)

Create a tensor grid for a rectangle from its 1d projections `x` and `y`.
Returns `unvec` to reshape a vector into a 2d array, `mtx` to evaluate a
function on the grid, `X`, `Y` to give the grid coordinates, and boolean array
`is_boundary` to identify the boundary points.
"""
function tensorgrid(x, y)
    m, n = length(x) - 1, length(y) - 1
    unvec(u) = reshape(u, m+1, n+1)
    mtx(h) = [h(x, y) for x in x, y in y]
    X = mtx((x, y) -> x)
    Y = mtx((x, y) -> y)
    is_boundary = trues(m+1, n+1)
    is_boundary[2:m, 2:n] .= false
    return mtx, X, Y, unvec, is_boundary
end
# end tensorgrid

# begin poissonfd
"""
    poissonfd(f, g, m, xspan, n, yspan)

Solve Poisson's equation on a rectangle by finite differences.
Function `f` is the forcing function and function `g` gives the
Dirichlet boundary condition. The rectangle is the tensor product of
intervals `xspan` and `yspan`,  and the discretization uses `m`+1
and `n`+1 points in the two coordinates.

Returns vectors defining the grid and a matrix of grid solution values.
"""
function poissonfd(f, g, m, xspan, n, yspan)
    # Discretize the domain.
    x, Dx, Dxx = FNC.diffmat2(m, xspan)
    y, Dy, Dyy = FNC.diffmat2(n, yspan)
    mtx, X, Y, unvec, is_boundary = tensorgrid(x, y)
    N = (m+1) * (n+1)   # total number of unknowns

    # Form the collocated PDE as a linear system.
    A = kron(I(n+1), sparse(Dxx)) + kron(sparse(Dyy), I(m+1))
    b = vec(mtx(f))

    # Apply Dirichlet condition.
    scale = maximum(abs, A[n+2, :])
    idx = vec(is_boundary)
    A[idx, :] = scale * I(N)[idx, :]        # Dirichet assignment
    b[idx] = scale * g.(X[idx], Y[idx])    # assigned values

    # Solve the linear system and reshape the output.
    u = A \ b
    return x, y, unvec(u)
end
# end poissonfd

# begin elliptic
"""
    elliptic(ϕ, g, m, xspan, n, yspan)

Solve the elliptic PDE
    `ϕ`(x, y, u, u_x, u_xx, u_y, u_yy) = 0
on the rectangle `xspan`x`yspan`, subject to `g`(x,y)=0 on the boundary. Uses
`m`+1 points in x by `n`+1 points in y in a Chebyshev discretization. Returns
vectors defining the grid and a matrix of grid solution values.
"""
function elliptic(ϕ, g, m, xspan, n, yspan)
    # Discretize the domain.
    x, Dx, Dxx = diffcheb(m, xspan)
    y, Dy, Dyy = diffcheb(n, yspan)
    mtx, X, Y, unvec, is_boundary = tensorgrid(x, y)
    N = (m+1) * (n+1)   # total number of unknowns

    # Identify boundary locations and evaluate the boundary condition.
    idx = vec(is_boundary)
    gb = g.(X[idx], Y[idx])

    # Evaluate the PDE+BC residual.
    function residual(u)
        U = unvec(u)
        R = ϕ(X, Y, U, Dx * U, Dxx * U, U * Dy', U * Dyy')    # PDE
        @. R[idx] = u[idx] - gb                               # boundary residual
        return vec(R)
    end

    # Solve the equation.
    u = levenberg(residual, vec(zeros(size(X))))[end]
    U = unvec(u)

    return function (ξ, η)
        v = [chebinterp(x, u, ξ) for u in eachcol(U)]
        return chebinterp(y, v, η)
    end
end
# end elliptic

"Evaluate Chebyshev interpolant with nodes x, values v, at point ξ"
function chebinterp(x, v, ξ)
    n = length(x) - 1
    w = (-1.0) .^ (0:n)
    w[[1, n+1]] .*= 0.5

    if ξ in x    # exactly at a node
        idx = findfirst(ξ .== x)
        return v[idx]
    else
        terms = @. w / (ξ - x)
        return sum(v .* terms) / sum(terms)
    end
end
