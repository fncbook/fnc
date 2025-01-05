---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-diffadv-wave)

We start with the discretization and initial condition.

```{code-cell}
m, n = 40, 40
x, Dx, Dxx = FNC.diffcheb(m, [-2, 2])
y, Dy, Dyy = FNC.diffcheb(n, [-2, 2])
mtx, X, Y, _ = FNC.tensorgrid(x, y)
U₀ = mtx( (x, y) -> (x + 0.2) * exp(-12 * (x^2 + y^2)) )
V₀ = zeros(size(U₀));
```

Note that because $u$ is known on the boundary, while $v$ is unknown over the full grid, there are two different sizes of vec/unvec operations. We also need to define functions to pack grid unknowns into a vector and to unpack them. When the unknowns for $u$ are packed, the boundary values are chopped off, and these are restored when unpacking.

```{code-cell}
_, _, _, unvec_v, _ = FNC.tensorgrid(x, y)
_, _, _, unvec_u, _ = FNC.tensorgrid(x[2:m], y[2:n])
chop = U -> U[2:m, 2:n]
extend = U -> [zeros(m+1) [zeros(1, n-1); U; zeros(1, n-1)] zeros(m+1)]
pack = (U, V) -> [vec(chop(U)); vec(V)]
N = (m-1) * (n-1)    # number of interior unknowns
unpack = w -> ( extend(unvec_u(w[1:N])), unvec_v(w[N+1:end]) )
```

We can now define and solve the IVP. Since this problem is hyperbolic, not parabolic, a nonstiff integrator is faster than a stiff one.

```{code-cell}
using OrdinaryDiffEq
function dw_dt(w, c, t)
    U, V = unpack(w)
    du_dt = V
    dv_dt = c^2 * (Dxx * U + U * Dyy')
    return pack(du_dt, dv_dt)
end

IVP = ODEProblem(dw_dt, pack(U₀, V₀), (0, 4.0), 1)
sol = solve(IVP, Tsit5())
U = t -> unpack(sol(t))[1]
```

```{code-cell}
anim = @animate for t in 0:4/100:4
    Ut = U(t)
    surface(x, y, Ut';
        layout=(1, 2), size=(640, 320),
        xlabel=L"x",  ylabel=L"y",  zaxis=((-0.1, 0.1), L"u(x,y)"),
        color=:redsblues,  alpha=0.66,  clims=(-0.1, 0.1), colorbar=:none,
        title="Wave equation",  dpi=150)
    contour!(x, y, Ut'; 
        levels=24,  subplot=2, 
        aspect_ratio=1,
        xlabel=L"x",  ylabel=L"y",
        color=:redsblues,  clims=(-0.1, 0.1), 
        colorbar=:none,  title=@sprintf("t = %.2f", t))
end
closeall();
mp4(anim, "diffadv-wave.mp4");
```

![Wave equation in 2d](diffadv-wave.mp4)
