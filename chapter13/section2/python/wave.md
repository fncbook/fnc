---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-diffadv-wave)

We start with the discretization and initial condition.

```{code-cell}
m, n = 40, 42
x, Dx, Dxx = FNC.diffcheb(m, [-2, 2])
y, Dy, Dyy = FNC.diffcheb(n, [-2, 2])
mtx, X, Y, _, _, _ = FNC.tensorgrid(x, y)

U0 = mtx(lambda x, y: (x + 0.2) * exp(-12 * (x**2 + y**2)))
V0 = zeros(U0.shape)
```

Note that because $u$ is known on the boundary, while $v$ is unknown over the full grid, there are two different sizes of vec/unvec operations. We also need to define functions to pack grid unknowns into a vector and to unpack them. When the unknowns for $u$ are packed, the boundary values are chopped off, and these are restored when unpacking.

```{code-cell}
_, _, _, vec_v, unvec_v, _ = FNC.tensorgrid(x, y)
_, _, _, vec_u, unvec_u, _ = FNC.tensorgrid(x[1:-1], y[1:-1])

def extend(U):
    UU = zeros((m+1, n+1))
    UU[1:-1, 1:-1] = U
    return UU

def chop(U):
    return U[1:-1, 1:-1]

def pack(U, V): 
    return hstack([vec_u(chop(U)), vec_v(V)])

N = (m-1) * (n-1)
def unpack(w):
    U = extend(unvec_u(w[:N]))
    V = unvec_v(w[N:])
    return U, V
```

We can now define and solve the IVP. Since this problem is hyperbolic, not parabolic, a nonstiff integrator is faster than a stiff one.

```{code-cell}
from scipy.integrate import solve_ivp
def dw_dt(t, w):
    U, V = unpack(w)
    dU_dt = V
    dV_dt = Dxx @ U + U @ Dyy.T
    return pack(dU_dt, dV_dt)

from scipy.integrate import solve_ivp
sol = solve_ivp(dw_dt, (0, 4), pack(U0, V0), method="RK45", dense_output=True)
U = lambda t: unpack(sol.sol(t))[0]
```

```{code-cell}
from matplotlib.animation import FuncAnimation
fig, ax = subplots()
obj = ax.pcolormesh(X, Y, U(0), vmin=-0.1, vmax=0.1, cmap="RdBu", shading="gouraud")
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel("$x$"),  ax.set_ylabel("$y$")
ax.set_aspect("equal")
ax.set_title("Wave equation in 2d")

def snapshot(t):
    global obj
    obj.remove()
    obj = ax.pcolormesh(X, Y, U(t), vmin=-0.1, vmax=0.1, cmap="RdBu", shading="gouraud")
    time_text.set_text(f"t = {t:.2f}")

anim = FuncAnimation(fig, snapshot, frames=linspace(0, 4, 91))
anim.save("wave-2d.mp4", fps=30);
close()
```

![Wave equation in 2d](wave-2d.mp4)
