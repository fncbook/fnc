---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
# Chapter 12

## Examples

```{code-cell}
:tags: [remove-output]
import Pkg; Pkg.activate("/Users/driscoll/Documents/GitHub/fnc")
using FundamentalsNumericalComputation
FNC.init_format()
```

### 12.1 @section-advection-traffic

(demo-traffic-advection-julia)=
``````{dropdown} @demo-traffic-advection

 In the following definition we allow the velocity $c$ to be specified as a parameter in the `ODEProblem`.

```{code-cell}
x, Dₓ, Dₓₓ = FNC.diffper(300, [-4, 4]);
f = (u, c, t) -> -c * (Dₓ*u);
```

The following initial condition isn't mathematically periodic, but the deviation is less than machine precision. We specify RK4 as the solver.  

```{code-cell}
u_init = @. 1 + exp( -3x^2 )
IVP = ODEProblem(f, u_init, (0., 4.), 2)
sol = solve(IVP, RK4());
```

```{code-cell}
:tags: hide-input
plt = plot(legend=:bottomleft,
    title="Advection with periodic boundary",
    xaxis=(L"x"),  yaxis=(L"u(x,t)")
    )
for t in (0:4) * 2/3
    plot!(x, sol(t), label=@sprintf("t=%.1f", t))
end
plt
```

An animation shows the solution nicely. The bump moves with speed 2 to the right, reentering on the left as it exits to the right because of the periodic conditions. 

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 4, 120) 
    plot(x, sol(t),
        title=@sprintf("Advection equation, t = %.2f", t),
        xaxis=(L"x"),  yaxis=([1, 2], L"u(x,t)"),
        dpi=150,
    )
end
mp4(anim,"figures/advection.mp4")
```

``````

(demo-traffic-solve-julia)=
``````{dropdown} @demo-traffic-solve
The following are parameters and a function relevant to defining the problem. 

```{code-cell}
ρc = 1080;  ρm = 380;  q_m = 10000;
dQ0 = ρ -> 4q_m * ρc^2 * (ρc-ρm) * ρm * (ρm-ρ) / (ρ*(ρc-2*ρm) + ρc*ρm)^3;
```

Here we create a discretization on $m=800$ points.

```{code-cell}
x, Dₓ, Dₓₓ = FNC.diffper(800, [0, 4]);
```

Next we define the ODE resulting from the method of lines.

```{code-cell}
ode = (ρ, ϵ, t) -> -dQ0.(ρ) .* (Dₓ*ρ) + ϵ * (Dₓₓ*ρ);
```

Our first initial condition has moderate density with a small bump. Because of the diffusion present, we use a stiff solver for the IVP.

```{code-cell}
ρ_init = @. 400 + 10 * exp( -20*(x-3)^2 )
IVP = ODEProblem(ode, ρ_init, (0., 1.), 0.02)
sol = solve(IVP, Rodas4P());
```

```{code-cell}
:tags: hide-input
plt = plot(legend=:topleft, 
    title="Traffic flow",
    xaxis=(L"x"),  yaxis=("car density"))
for t in 0:0.2:1
    plot!(x, sol(t), label=@sprintf("t=%.1f", t))
end
plt
```

The bump slowly moves backward on the roadway, spreading out and gradually fading away due to the presence of diffusion.

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0,0.9,91) 
    plot(x, sol(t),
        xaxis=(L"x"), yaxis=([400,410], "density"),
        dpi=150,    
        title=@sprintf("Traffic flow, t=%.2f",t) 
        )
end
mp4(anim,"figures/traffic-fade.mp4")
```

Now we use an initial condition with a larger bump. Note that the scale on the $y$-axis is much different for this solution.

```{code-cell}
ρ_init = @. 400 + 80 * exp( -16*(x - 3)^2 )
IVP = ODEProblem(ode, ρ_init, (0., 0.5), 0.02)
sol = solve(IVP, Rodas4P());
```

```{code-cell}
:tags: hide-input
plt = plot(legend=:topleft,
    title="Traffic jam",
    xaxis=(L"x"),  yaxis=("car density")
    )
for t in range(0, 5, 11)
    plot!(x, sol(t), label=@sprintf("t=%.1f", t))
end
plt
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 0.5, 101) 
    plot(x, sol(t);
        xaxis=(L"x"),  yaxis=([400,480], "density"),
        dpi=150,    
        title=@sprintf("Traffic jam, t=%.2f",t) 
        )
end
mp4(anim,"figures/traffic-jam.mp4")
```

In this case the density bump travels backward along the road. It also steepens on the side facing the incoming traffic and decreases much more slowly on the other side. A motorist would experience this as an abrupt increase in density, followed by a much more gradual decrease in density and resulting gradual increase in speed. (You also see some transient, high-frequency oscillations. These are caused by instabilities, as we discuss in simpler situations later in this chapter.)

``````

### 12.2 @section-advection-upwind
(demo-upwind-cfl-julia)=
``````{dropdown} @demo-upwind-cfl
For time stepping, we use the adaptive explicit method `RK4`.

```{code-cell}
function demo(m)
  x, Dₓ = FNC.diffper(m, [0, 1])
  uinit = @. exp( -80 * (x - 0.5)^2 )
  ode = (u, c, t) -> -c * (Dₓ*u)
  IVP = ODEProblem(ode, uinit, (0., 2.), 2.)
  return x, solve(IVP, RK4())
end
x,u = demo(400);
```

```{code-cell}
:tags: hide-input
t = range(0, 2, 81);
U = reduce(hcat, u(t) for t in t)
contour(x, t, U'; 
    color=:redsblues, 
    clims=(-1, 1),
    xaxis=(L"x"),  yaxis=(L"t"),
    title="Linear advection",
    right_margin=3Plots.mm
    )
```

In the space-time plot above, you can see the initial hump traveling rightward at constant speed. It fully traverses the domain once for each integer multiple of $t=1/2$. 

If we cut $h$ by a factor of 2 (i.e., double $m$), then the CFL condition suggests that the time step should be cut by a factor of 2 also.

```{code-cell}
println("Number of time steps for m = 400: $(length(u.t))")
x, u = demo(800)
println("Number of time steps for m = 800: $(length(u.t))")
```
``````

(demo-upwind-direction-julia)=
``````{dropdown} @demo-upwind-direction
If we solve advection over $[0,1]$ with velocity $c=-1$, the right boundary is in the upwind/inflow direction. Thus a well-posed boundary condition is $u(1,t)=0$.

We'll pattern a solution after {numref}`Function {number} <function-parabolic>`. Since $u(x_m,t)=0$, we define the ODE interior problem {eq}`mol-interior` for $\mathbf{v}$ without $u_m$. For each evaluation of $\mathbf{v}'$, we must extend the data back to $x_m$ first.

```{code-cell}
m = 100
x, Dₓ = FNC.diffmat2(m, [0, 1])

interior = 1:m
extend = v -> [v; 0]

function ode!(f, v, c, t)
    u = extend(v)
    uₓ = Dₓ * u
    @. f = -c * uₓ[interior]
end;
```

Now we solve for an initial condition that has a single hump.

```{code-cell}
init = @. exp( -80*(x[interior] - 0.5)^2 )
ivp = ODEProblem(ode!, init, (0., 1), -1)
u = solve(ivp);
```

```{code-cell}
t = range(0, 0.75, 80)
U = reduce(hcat, extend(u(t)) for t in t)
contour(x, t, U';
    color=:blues,  clims=(0, 1), 
    xaxis=(L"x"),  yaxis=(L"t"),
    title="Advection with inflow BC"
    )
```

We find that the hump gracefully exits out the downwind end.

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 1, 161) 
    plot(x, extend(u(t));
        label=@sprintf("t = %.4f", t), 
        xaxis=(L"x"),  yaxis=(L"u(x, t)", (0, 1)), 
        title="Advection equation with inflow BC",
        dpi=150
        )
end
mp4(anim,"figures/upwind-inflow.mp4")
```

If instead of $u(1,t)=0$ we were to try to impose the downwind condition $u(0,t)=0$, we only need to change the index of the interior nodes and where to append the zero value.

```{code-cell}
interior = 2:m+1
extend = v -> [0; v]

init = @. exp( -80*(x[interior] - 0.5)^2 )
ivp = ODEProblem(ode!, init, (0., 0.25), -1)
u = solve(ivp);
```

```{code-cell}
:tags: hide-input
t = range(0, 0.2, 61)
U = reduce(hcat, extend(u(t)) for t in t)
contour(x, t, U'; 
    color=:redsblues,  clims=(-1, 1),
    xaxis=(L"x"),  yaxis=(L"t"), 
    title="Advection with outflow BC",
    right_margin=3Plots.mm)
```

This time, the solution blows up as soon as the hump runs into the boundary because there are conflicting demands there.

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0, 0.2, 41) 
    plot(x, extend(u(t));
        label=@sprintf("t = %.4f", t),
        xaxis=(L"x"),  yaxis=(L"u(x,t)", (0, 1)), 
        title="Advection equation with outflow BC",dpi=150)
end
mp4(anim,"figures/upwind-outflow.mp4")
```
``````

### 12.3 @section-advection-absstab

(demo-absstab-advection-julia)=
``````{dropdown} @demo-absstab-advection
For $c=1$ we get purely imaginary eigenvalues.

```{code-cell}
:tags: hide-input
x,Dₓ = FNC.diffper(40, [0, 1])
λ = eigvals(Dₓ);

scatter(real(λ), imag(λ);
    aspect_ratio = 1, 
    xlabel="Re λ",  ylabel="Im λ", 
    frame=:zerolines, 
    title="Eigenvalues for pure advection", 
    leg=:none
    )
```

Let's choose a time step of $\tau=0.1$ and compare to the stability regions of the Euler and backward Euler time steppers (shown as shaded regions):

```{code-cell}
:tags: hide-input
zc = @. cispi(2 * (0:360) / 360);     # points on |z|=1
z = zc .- 1;                          # shift left by 1
plot(Shape(real(z), imag(z)), color=RGB(.8, .8, 1))

ζ = 0.1 * λ
scatter!(real(ζ), imag(ζ);
    aspect_ratio=1,
    xaxis=("Re ζ", [-5, 5]),
    yaxis=("Im ζ", [-5, 5]),
    frame=:zerolines,
    title="Euler for advection"
    )
```

In the Euler case it's clear that *no* real value of $\tau>0$ is going to make $\zeta$ values fit within the stability region. Any method whose stability region includes none of the imaginary axis is an unsuitable choice for advection.

```{code-cell}
:tags: hide-input
z = zc .+ 1;                        # shift circle right by 1
plot(Shape([-6, 6, 6, -6], [-6, -6, 6, 6]), color=RGB(.8, .8, 1))
plot!(Shape(real(z), imag(z)), color=:white)

scatter!(real(ζ), imag(ζ);
    aspect_ratio=1,
    xaxis=("Re ζ", [-5, 5]),
    yaxis=("Im ζ", [-5, 5]),
    frame=:zerolines,
    title="Backward Euler for advection"
    )
```

The A-stable backward Euler time stepping tells the exact opposite story: it will be absolutely stable for any choice of the time step $\tau$.
``````

(demo-absstab-advdiff-julia)=
``````{dropdown} @demo-absstab-advdiff
The eigenvalues of advection-diffusion are near-imaginary for $\epsilon\approx 0$ and get closer to the negative real axis as $\epsilon$ increases.

```{code-cell}
:tags: hide-input
plt = plot(leg=:topleft,
    aspect_ratio=1,
    xlabel="Re ζ",  ylabel="Im ζ",
    title="Eigenvalues for advection-diffusion"
    )
x, Dₓ, Dₓₓ = FNC.diffper(40, [0, 1]);
for ϵ in [0.001, 0.01, 0.05]
    λ = eigvals(-Dₓ + ϵ*Dₓₓ)
    scatter!(real(λ), imag(λ), m=:o, label="\\epsilon = $ϵ")
end
plt
```
``````

(demo-absstab-inflow-julia)=
``````{dropdown} @demo-absstab-inflow
Deleting the last row and column places all the eigenvalues of the discretization into the left half of the complex plane. 

```{code-cell}
x, Dₓ, _ = FNC.diffcheb(40, [0, 1])
A = Dₓ[1:end-1, 1:end-1];     # delete last row and column
λ = eigvals(A);
```

```{code-cell}
:tags: hide-input

scatter(real(λ), imag(λ);
    m=3,  aspect_ratio=1,
    label="", 
    xaxis=([-300, 100], "Re λ"), 
    yaxis=("Im λ"),
    title="Eigenvalues of advection with zero inflow") 
```

Note that the rightmost eigenvalues have real part at most

```{code-cell}
maximum(real(λ))
```

Consequently all solutions decay exponentially to zero as $t\to\infty$. This matches our observation of the solution: eventually, everything flows out of the domain.

``````

### 12.4 @section-advection-wave

(demo-wave-boundaries-julia)=
``````{dropdown} @demo-wave-boundaries

```{code-cell}
m = 200
x, Dₓ = FNC.diffcheb(m, [-1, 1]);
```

The boundary values of $u$ are given to be zero, so they are not unknowns in the ODEs. Instead they are added or removed as necessary.

```{code-cell}
extend = v -> [0; v; 0]
chop = u -> u[2:m];
```

The following function computes the time derivative of the system at interior points.

```{code-cell}
ode = function(w, c, t)
    u = extend(w[1:m-1])
    z = w[m:2m]
    dudt = Dₓ * z
    dzdt = c^2 * (Dₓ * u)
    return [ chop(dudt); dzdt ]
end;
```

Our initial condition is a single hump for $u$.

```{code-cell}
u_init = @. exp( -100*(x + 0.5)^2 )
z_init = -u_init
w_init = [ chop(u_init); z_init ];  
```

Because the wave equation is hyperbolic, we can use a nonstiff explicit solver.

```{code-cell}
IVP = ODEProblem(ode, w_init ,(0., 2.), 2)
w = solve(IVP, RK4());
```

We plot the results for the original $u$ variable only. Its interior values are at indices `1:m-1` of the composite $\mathbf{w}$ variable.

```{code-cell}
t = range(0, 2, 80)
U = [extend(w(t)[1:m-1]) for t in t]
contour(x, t, hcat(U...)';
    color=:redsblues,  clims=(-1, 1),
    levels=24,
    xlabel=L"x",  ylabel=L"t",
    title="Wave equation",
    right_margin=3Plots.mm
    )
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0 ,2, 120)
    plot(x, extend(w(t)[1:m-1]);
        label=@sprintf("t=%.3f",t),
        xaxis=(L"x"),  yaxis=([-1, 1], L"u(x,t)"),
        dpi=150,    
        title="Wave equation"
        )
end
mp4(anim,"figures/wave-boundaries.mp4")
```

The original hump breaks into two pieces of different amplitudes, each traveling with speed $c=2$. They pass through one another without interference. When a hump encounters a boundary, it is perfectly reflected, but with inverted shape. At time $t=2$, the solution looks just like the initial condition.

``````

(demo-wave-speed-julia)=
``````{dropdown} @demo-wave-speed
The ODE implementation has to change slightly.

```{code-cell}
ode = function(w,c,t)
    u = extend(w[1:m-1])
    z = w[m:2m]
    dudt = Dₓ*z
    dzdt = c.^2 .* (Dₓ*u)
    return [ chop(dudt); dzdt ]
end;
```

The variable wave speed is passed as an extra parameter through the IVP solver.

```{code-cell}
c = @. 1 + (sign(x)+1)/2
IVP = ODEProblem(ode, w_init, (0., 5.), c)
w = solve(IVP, RK4());
```

```{code-cell}
t = range(0, 5, 80)
U = [extend(w(t)[1:m-1]) for t in t]
contour(x, t, hcat(U...)';
    color=:redsblues,  clims=(-1,1),
    levels=24,
    xlabel=L"x",  ylabel=L"t",
    title="Wave equation",
    right_margin=3Plots.mm
    )
```

```{code-cell}
:tags: hide-input
anim = @animate for t in range(0,5,181)
    plot(Shape([-1, 0, 0, -1], [-1, -1, 1, 1]), color=RGB(.8, .8, .8), l=0, label="")
    plot!(x, extend(w(t, idxs=1:m-1));
        label=@sprintf("t=%.2f", t), 
        xaxis=(L"x"),  yaxis=([-1, 1], L"u(x,t)"),
        dpi=150,   
        title="Wave equation, variable speed" 
        )
end
mp4(anim,"figures/wave-speed.mp4")
```

Each pass through the interface at $x=0$ generates a reflected and transmitted wave. By conservation of energy, these are both smaller in amplitude than the incoming bump.
``````