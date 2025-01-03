---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 13

## Functions

(function-tensorgrid-matlab)=
``````{dropdown} Create a tensor-product grid
```{literalinclude} fnc/tensorgrid.m
:language: matlab
:linenos: true
```
``````

(function-poissonfd-matlab)=
``````{dropdown} Solution of Poisson's equation by finite differences
```{literalinclude} fnc/poissonfd.m
:language: matlab
:linenos: true
```
``````

(function-elliptic-matlab)=
``````{dropdown} Solution of elliptic PDE by Chebyshev collocation
```{literalinclude} fnc/elliptic.m
:language: matlab
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: remove-cell
addpath /Users/driscoll/Documents/GitHub/fnc/matlab/fnc
addpath /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```

### 13.1 @section-twodim-tensorprod

(demo-tensorprod-gridfun-matlab)=
``````{dropdown} @demo-tensorprod-gridfun
Here is the grid from {numref}`Example {number} <example-tensorprod-smallgrid>`.

```{code-cell}
m = 4;   
x = linspace(0, 2, m+1);
n = 2;   
y = linspace(1, 3, n+1);
```

For a given $f(x,y)$ we can find $\operatorname{mtx}(f)$ by using a comprehension syntax.

```{code-cell}
[mtx, X, Y] = tensorgrid(x, y);
f = @(x, y) cos(pi * x .* y - y);
F = mtx(f)
```

::::{grid} 1 1 2 2
We can make a nice plot of the function by first choosing a much finer grid. However, the contour and surface plotting functions expect the *transpose* of mtx($f$).
:::{card}
To emphasize departures from a zero level, use a colormap such as `redsblues` and set the color limits to be symmetric around zero.
:::
::::

::::{warning}
The contour and surface plotting functions expect the *transpose* of the outputs of `mtx`. If you forget to do that, the $x$ and $y$ axes will be swapped.
::::

```{code-cell}
m = 80;  x = linspace(0, 2, m+1);
n = 60;  y = linspace(1, 3, n+1);
[mtx, X, Y] = tensorgrid(x, y);
F = mtx(f);

pcolor(X', Y', F')
shading interp
colormap(redsblues),  colorbar
axis equal
xlabel("x"),  ylabel("y")
```
``````

(demo-tensorprod-disksphere-matlab)=
``````{dropdown} @demo-tensorprod-disksphere
For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = linspace(0, 1, 41);
theta = linspace(0, 2*pi, 121);
[mtx, R, Theta] = tensorgrid(r, theta);
F = mtx(@(r, theta) 1 - r.^4);
clf,  colormap(parula)
contourf(R', Theta', F', 20)
shading flat
xlabel("r"),  ylabel("\theta"), 
title("A polar function")   
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. For this we create matrices for the coordinate functions $x$ and $y$.

```{code-cell}
X = R .* cos(Theta);  Y = R .* sin(Theta);
contourf(X', Y', F', 20)
axis equal,  shading interp  
xlabel('x'),  ylabel('y')
title('Function over the unit disk')  
```

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

``````

(demo-tensorprod-diff-matlab)=
``````{dropdown} @demo-tensorprod-diff
We define a function and, for reference, its two exact partial derivatives.

```{code-cell}
u = @(x, y) sin(pi * x .* y - y);
du_dx = @(x, y) pi * y .* cos(pi * x .* y - y);
du_dy = @(x, y) (pi * x - 1) .* cos(pi * x .* y - y);
```

We will use an equispaced grid and second-order finite differences as implemented by `diffmat2`. First, we have a look at a plots of the exact partial derivatives.

```{code-cell}
m = 80;  [x, Dx] = diffmat2(m, [0, 2]);
n = 60;  [y, Dy] = diffmat2(n, [1, 4]);
[mtx, X, Y] = tensorgrid(x, y);
U = mtx(u);
dU_dX = mtx(du_dx);
dU_dY = mtx(du_dy);
```

```{code-cell}
clf,  subplot(1, 2, 1)
pcolor(X', Y', dU_dX')
axis equal,  shading interp
title('∂u/∂x')
subplot(1, 2, 2)
pcolor(X', Y', dU_dY')
axis equal,  shading interp
title('∂u/∂y')
```

Now we compare the exact partial derivatives with their finite-difference approximations. Since these are signed errors, we use a colormap that is symmetric around zero.

```{code-cell}
err = dU_dX - Dx * U;
subplot(1, 2, 1)
pcolor(X', Y', err')
colorbar,  clim([-.4, .4])
axis equal,  shading interp
title('error in ∂u/∂x')

err = dU_dY - U * Dy';
subplot(1,2,2)
pcolor(X', Y', err')
colorbar,  clim([-.1, .1])
axis equal,  shading interp
colormap(redsblues)
title('error in ∂u/∂y')
```

Not surprisingly, the errors are largest where the derivatives themselves are largest.
``````

### 13.2 @section-twodim-diffadv

(demo-diffadv-vec-matlab)=
``````{dropdown} @demo-diffadv-vec

```{code-cell}
m = 4;  n = 3;
x = linspace(0, 2, m+1);
y = linspace(-3, 0, n+1);

f = @(x, y) cos(0.75*pi * x .* y - 0.5*pi * y);
[mtx, X, Y, vec, unvec] = tensorgrid(x, y);
F = mtx(f);
disp("function on a 4x3 grid:")
disp(F)
```

```{code-cell}
disp("vec(F):")
disp(vec(F))
```

The `unvec` operation is the inverse of vec.

```{code-cell}
disp("unvec(vec(F)):")
disp(unvec(vec(F)))
```
``````

(demo-diffadv-heat-matlab)=
``````{dropdown} @demo-diffadv-heat

```{code-cell}
m = 60;  n = 40;
[x, Dx, Dxx] = diffper(m, [-1, 1]);
[y, Dy, Dyy] = diffper(n, [-1, 1]);
[mtx, X, Y, vec, unvec] = tensorgrid(x, y);
```

Note that the initial condition should also be periodic on the domain.

```{code-cell}
U0 = sin(4*pi*X) .* exp( cos(pi*Y) );
clf,  surf(X', Y', U0')
mx = max(abs(vec(U0)));
clim([-mx, mx]),  shading interp
colormap(redsblues)
xlabel('x'),  ylabel('y')  
title('Initial condition')  
```

This function computes the time derivative for the unknowns. The actual calculations take place using the matrix shape.

```{literalinclude} f13_2_heat.m
:language: matlab
```

Since this problem is parabolic, a stiff integrator is appropriate.

```{code-cell}
ivp = ode(ODEFcn=@f13_2_heat);
ivp.InitialTime = 0;
ivp.InitialValue = vec(U0);
ivp.Parameters = {0.1, Dxx, Dyy, vec, unvec};
ivp.Solver = "stiff";
sol = solutionFcn(ivp, 0, 0.2);
U = @(t) unvec(sol(t));
```

```{code-cell}
:tags: hide-input
surf(X', Y', U0')
mx = max(abs(vec(U0)));
clim([-mx, mx]),  shading interp
colormap(redsblues)
xlabel('x'),  ylabel('y')  
title('Initial condition') 
```

::::{grid} 1 1 2 2
Here is an animation of the solution.
:::{card}
Here `clims` are set so that colors remain at fixed values throughout the animation.
:::
::::

```{code-cell}
title('Heat equation on a periodic domain')
vid = VideoWriter("figures/2d-heat.mp4","MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 0.2, 61)
    cla, surf(X', Y', U(t)')
    shading interp
    str = sprintf("t = %.2f", t);
    text(-0.9, 0.75, 2, str, fontsize=14);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

![Heat equation in 2d](figures/2d-heat.mp4)

``````

(demo-diffadv-advdiff-matlab)=
``````{dropdown} @demo-diffadv-advdiff

The first step is to define a discretization of the domain and the initial state.

```{code-cell}
m = 50;  n = 40;
[x, Dx, Dxx] = diffcheb(m, [-1, 1]);
[y, Dy, Dyy] = diffcheb(n, [-1, 1]);
[mtx, X, Y] = tensorgrid(x, y);
u_init = @(x, y) (1+y) .* (1-x).^4 .* (1+x).^2 .* (1-y.^4);
```

There are really two grids now: the full grid and the subset grid of interior points. Since the IVP unknowns are on the interior grid, that is the one we need to change shapes on. We also need the functions `extend` and `chop` to add and remove boundary values.

```{code-cell}
[~, ~, ~, vec, unvec] = tensorgrid(x(2:m), y(2:n));
chop = @(U) U(2:m, 2:n);
z = zeros(1, n-1);
extend = @(U) [ zeros(m+1, 1) [z; U; z] zeros(m+1, 1)];
pack = @(U) vec(chop(U));
unpack = @(u) extend(unvec(u));
```

Now we can define and solve the IVP using a stiff solver.

```{literalinclude} f13_2_advdiff.m
:language: matlab
```

```{code-cell}
ivp = ode(ODEFcn=@f13_2_advdiff);
ivp.InitialTime = 0;
ivp.InitialValue = pack(mtx(u_init));
ivp.Parameters = {0.05, Dx, Dxx, Dy, Dyy, pack, unpack};
ivp.Solver = "stiff";
sol = solutionFcn(ivp, 0, 2);
```

When we evaluate the solution at a particular value of $t$, we get a vector of the interior grid values. The same `unpack` function above converts this to a complete matrix of grid values.

```{code-cell}
U = @(t) unpack(sol(t));

clf,  pcolor(X', Y', U(0.5)')
clim([0, 2]), shading interp
axis equal,  colormap(sky), colorbar
title('Advection-diffusion at t = 0.5')  
xlabel('x'),  ylabel('y')  
```

```{code-cell}
:tags: hide-input
hold on
vid = VideoWriter("figures/2d-advdiff.mp4","MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 2, 81)
    cla, pcolor(X', Y', U(t)')
    shading interp
    str = sprintf("t = %.2f", t);
    text(-1.5, 0.75, str, fontsize=14);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

![Advection-diffusion in 2d](figures/2d-advdiff.mp4)
``````

(demo-diffadv-wave-matlab)=
``````{dropdown} @demo-diffadv-wave
We start with the discretization and initial condition.

```{code-cell}
m = 40; n = 42;
[x, Dx, Dxx] = diffcheb(m, [-2, 2]);
[y, Dy, Dyy] = diffcheb(n, [-2, 2]);
[mtx, X, Y] = tensorgrid(x, y);

u_init = @(x, y) (x+0.2) .* exp(-12*(x.^2 + y.^2));
U0 = mtx(u_init);
V0 = zeros(size(U0));
```

Note that because $u$ is known on the boundary, while $v$ is unknown over the full grid, there are two different sizes of vec/unvec operations. We also need to define functions to pack grid unknowns into a vector and to unpack them. When the unknowns for $u$ are packed, the boundary values are chopped off, and these are restored when unpacking.

```{code-cell}
[~, ~, ~, vec_v, unvec_v] = tensorgrid(x, y);
[~, ~, ~, vec_u, unvec_u] = tensorgrid(x(2:m), y(2:n));

chop = @(U) U(2:m, 2:n);
z = zeros(1, n-1);
extend = @(U) [ zeros(m+1, 1) [z; U; z] zeros(m+1, 1)];
pack = @(U, V) [vec_u(chop(U)); vec_v(V)];
N = (m-1) * (n-1);
unpack = @(u) f13_2_wave_unpack(u, N, unvec_u, unvec_v, extend);
```

```{literalinclude} f13_2_wave_unpack.m
:language: matlab
```


We can now define and solve the IVP. Since this problem is hyperbolic, not parabolic, a nonstiff integrator is faster than a stiff one.

```{literalinclude} f13_2_wave.m
:language: matlab
```

```{code-cell}
ivp = ode(ODEFcn=@f13_2_wave);
ivp.InitialTime = 0;
ivp.InitialValue = pack(U0, V0);
ivp.Parameters = {Dxx, Dyy, pack, unpack};
ivp.Solver = "nonstiff";
sol = solutionFcn(ivp, 0, 4);
```

```{code-cell}
clf
[U, V] = unpack(sol(0.5));
pcolor(X', Y', U')
axis equal,  clim([-0.1, 0.1])
colormap(redsblues),  shading interp
xlabel("x"),  ylabel("y")
title("Wave equation at t = 0.5")
```

```{code-cell}
:tags: hide-input
hold on
vid = VideoWriter("figures/2d-wave.mp4","MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 4, 121)
    [U, V] = unpack(sol(t));
    cla, pcolor(X, Y, U)
    shading interp
    str = sprintf("t = %.2f", t);
    text(-3, 1.75, str, fontsize=14);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
  
```

![Wave equation in 2d](figures/2d-wave.mp4)
``````

### 13.3 @section-twodim-laplace

(demo-laplace-kron-matlab)=
``````{dropdown} @demo-laplace-kron

```{code-cell}
A = [1, 2; -2, 0];
B = [1, 10, 100; -5, 5, 3];
disp("A:")
disp(A)
disp("B:")
disp(B)
```

Applying the definition manually, we get

```{code-cell}
A_kron_B = [
    A(1,1)*B  A(1,2)*B;
    A(2,1)*B  A(2,2)*B
    ]
```

```{index} ! MATLAB; kron
```

But it makes more sense to use `kron`.

```{code-cell}
kron(A, B)
```
``````

(demo-laplace-fd-matlab)=
``````{dropdown} @demo-laplace-fd
We make a crude discretization for illustrative purposes.

```{code-cell}
m = 6;  n = 5;
[x, Dx, Dxx] = diffmat2(m, [0, 3]);
[y, Dy, Dyy] = diffmat2(n, [-1, 1]);
[mtx, X, Y, vec, unvec, is_boundary] = tensorgrid(x, y);
```

Next, we define $\phi$ and evaluate it on the grid to get the forcing vector of the linear system.

```{code-cell}
phi = @(x, y) x.^2 - y + 2;
b = vec(mtx(phi));
```

Here are the coefficients for the PDE collocation, before any modifications are made for the boundary conditions. The combination of Kronecker products and finite differences produces a characteristic sparsity pattern.

```{code-cell}
A = kron(speye(n+1), sparse(Dxx)) + kron(sparse(Dyy), speye(m+1));
clf,  spy(A)
title("System matrix before boundary conditions")
```

The number of equations is equal to $(m+1)(n+1)$, which is the total number of points on the grid.

```{code-cell}
N = length(b)
```

We now use the Boolean array that indicates where the boundary points lie in the grid.

```{code-cell}
spy(is_boundary)
title("Boundary points")
```

::::{grid} 1 1 2 2
In order to impose Dirichlet boundary conditions, we replace the boundary rows of the system by rows of the identity.
:::{card}
Changing rows of a sparse array requires that the operands be in a particular sparse representation called `lil`. The conversion isn't done automatically because it can be slow and you are encouraged to avoid it when possible. We're just trying to keep things conceptually simple here.
:::
::::

```{code-cell}
I = speye(N);
idx = vec(is_boundary);
A(idx, :) = I(idx, :);

spy(A)
title("System matrix with boundary conditions")
```

Finally, we must replace the rows in the vector $\mathbf{b}$ by the boundary values being assigned to the boundary points. Here, we let the boundary values be zero everywhere.

```{code-cell}
b(idx) = 0;
```

Now we can solve for $\mathbf{u}$ and reinterpret it as the matrix-shaped $\mathbf{U}$, the solution on our grid.

```{code-cell}
u = A \ b;
U = unvec(u)
```
``````

(demo-laplace-poisson-matlab)=
``````{dropdown} @demo-laplace-poisson

First we define the problem on $[0,1]\times[0,2]$.

```{code-cell}
f = @(x, y) -sin(3*x .* y - 4*y) .* (9*y.^2 + (3*x - 4).^2);
g = @(x, y) sin(3*x .* y - 4*y);
xspan = [0, 1];
yspan = [0, 2];
```

Here is the finite-difference solution.

```{code-cell}
[X, Y, U] = poissonfd(f, g, 40, xspan, 60, yspan);

clf, surf(X', Y', U')
colormap(parula),  shading interp
colorbar
title("Solution of Poisson's equation")
xlabel("x"),  ylabel("y"),  zlabel("u(x,y)")
```

Since this is an artificial problem with a known solution, we can plot the error, which is a smooth function of $x$ and $y$. It must be zero on the boundary; otherwise, we have implemented boundary conditions incorrectly.

```{code-cell}
err = g(X, Y) - U;
mx = max(abs(vec(err)));
pcolor(X', Y', err')
colormap(redsblues),  shading interp
clim([-mx, mx]),  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Error")
```
``````

### 13.4 @section-twodim-nonlinear

(demo-nonlinear2d-mems-matlab)=
``````{dropdown} @demo-nonlinear2d-mems
All we need to define are $\phi$ from {eq}`nonlinpdepde` for the PDE, and a trivial zero function for the boundary condition.

```{code-cell}
lambda = 1.5;
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) Uxx + Uyy - lambda ./ (U + 1).^2;
g = @(x, y) zeros(size(x));
```

Here is the solution for $m=15$, $n=8$.

```{code-cell}
u = elliptic(phi, g, 15, [0, 2.5], 8, [0, 1]);
disp(sprintf("solution at (2, 0.6) is %.7f", u(2, 0.6)))
```

```{code-cell}
:tags: hide-input
x = linspace(0, 2.5, 91);
y = linspace(0, 1, 51);
[mtx, X, Y] = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(flipud(sky)),  shading interp,  colorbar
axis equal
xlabel("x"),  ylabel("y")
title("Deflection of a MEMS membrane")
```

In the absence of an exact solution, how can we be confident that the solution is accurate? First, the Levenberg iteration converged without issuing a warning, so we should feel confident that the discrete equations were solved. Assuming that we encoded the PDE correctly, the remaining source of error is truncation from the discretization. We can estimate that by refining the grid a bit and seeing how much the numerical solution changes.

```{code-cell}
x_test = linspace(0, 2.5, 6);
y_test = linspace(0, 1 , 5);
mtx_test = tensorgrid(x_test, y_test);
format long
mtx_test(u)
```

```{code-cell}
u = elliptic(phi, g, 25, [0, 2.5], 14, [0, 1]);
mtx_test(u)
```

```{code-cell}
:tags: remove-cell
format
```

The original solution seems to be accurate to about four digits.

``````

(demo-nonlinear-advdiff-matlab)=
``````{dropdown} @demo-nonlinear-advdiff

```{code-cell}
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) 1 - Ux - 2*Uy + 0.05*(Uxx + Uyy);
g = @(x, y) zeros(size(x));
u = elliptic(phi, g, 32, [-1, 1], 32, [-1, 1]);
```

```{code-cell}
:tags: hide-input
x = linspace(-1, 1, 80);
y = x;
mtx = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(parula),  shading interp,  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Steady advection–diffusion")
```

``````

(demo-nonlinear2d-allencahn-matlab)=
``````{dropdown} @demo-nonlinear2d-allencahn

The following defines the PDE and a nontrivial Dirichlet boundary condition for the square $[0,1]^2$.

```{code-cell}
phi = @(X, Y, U, Ux, Uxx, Uy, Uyy) U .* (1-U.^2) + 0.05*(Uxx + Uyy);
g = @(x, y) tanh(5*(x + 2*y - 1));
```

We solve the PDE and then plot the result.

```{code-cell}
u = elliptic(phi, g, 36, [0, 1], 36, [0, 1]);
```

```{code-cell}
:tags: hide-input
x = linspace(0, 1, 80);
y = x;
mtx = tensorgrid(x, y);
clf,  pcolor(x, y, mtx(u)')
colormap(parula),  shading interp,  colorbar
axis equal,  xlabel("x"),  ylabel("y")
title("Steady Allen–Cahn")
```
``````

