---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 11

## Functions

(function-diffper-matlab)=
``````{dropdown} Differentiation matrices for periodic end conditions
:open:
```{literalinclude} FNC-matlab/diffper.m
:language:matlab
:linenos: true
```
``````

(function-parabolic-matlab)=
``````{dropdown} Solution of parabolic PDEs by the method of lines
:open:
```{literalinclude} FNC-matlab/parabolic.m
:language:matlab
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init; pwd;
```

### 11.1 @section-diffusion-blackscholes

(demo-blackscholes-solve-matlab)=
``````{dropdown} @demo-blackscholes-solve
:open:
We consider the Black–Scholes problem for the following parameter values:

```{code-cell}
Smax = 8;  T = 6;
K = 3;  sigma = 0.06;  r = 0.08;
```

We discretize space and time.

```{code-cell}
m = 200;  h = Smax / m;
x = h * (0:m)';
n = 1000;  tau = T / n;
t = tau * (0:n)';
lambda = tau / h^2;  mu = tau / h;
```

We set the initial condition and then march forward in time.

```{code-cell}
V = zeros(m+1, n+1);
V(:, 1) = max(0, x-K);
for j = 1:n
    % Fictitious value from Neumann condition.
    Vfict = 2*h + V(m, j);
    Vj = [ V(:, j); Vfict ];
    % First row is zero by the Dirichlet condition.
    for i = 2:m+1 
        diff1 = (Vj(i+1) - Vj(i-1));
        diff2 = (Vj(i+1) - 2*Vj(i) + Vj(i-1));
        V(i, j+1) = Vj(i) ...
            + (lambda * sigma^2* x(i)^2/2) * diff2  ...
            + (r*mu * x(i))/2 * diff1 - r*tau * Vj(i);
    end 
end
```

Here is a plot of the solution after every 250 time steps.

```{code-cell}
index_times = 1:250:n+1;
show_times = t(index_times);
clf
for j = index_times
    str = sprintf("t = %.2f", t(j));
    plot(x, V(:, j), displayname=str) 
    hold on
end
title('Black-Scholes solution')
xlabel('stock price'),  ylabel('option value')
axis tight,  grid on
legend(location="northwest")
```

```{index} MATLAB; animation
```

Alternatively, here is an animation of the solution.

```{code-cell}
:tag: remove-output
clf
plot(x, V(:,1))
hold on,  grid on
axis([0, 8, 0, 6])
title('Black-Scholes solution') 
xlabel('stock price'),  ylabel('option value')
vid = VideoWriter("figures/black-scholes-6.mp4","MPEG-4");
vid.Quality = 85;
open(vid)
for frame = 1:10:n+1
    cla, plot(x, V(:, frame))
    str = sprintf("t = %.2f", t(frame));
    text(0.4, 5.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

![Black–Scholes solution](figures/black-scholes-6.mp4)

The results are easy to interpret, recalling that the time variable really means *time until strike*. Say you are close to the option's strike time. If the current stock price is, say, $S=2$, then it's not likely that the stock will end up over the strike price $K=3$, and therefore the option has little value. On the other hand, if presently $S=3$, then there are good odds that the option will be exercised at the strike time, and you will need to pay a substantial portion of the stock price in order to take advantage. As the time to strike increases, there is an expectation that the stock price is more likely to rise somewhat, making the value of the option larger at each fixed $S$. 
``````

(demo-blackscholes-unstable-matlab)=
``````{dropdown} @demo-blackscholes-unstable
:open:
Let's try to do everything the same as in @demo-blackscholes-solve, but extending the simulation time to $T=8$.

```{code-cell}
T = 8;
n = 1000;  tau = T / n;
t = tau*(0:n)';
lambda = tau / h^2;  mu = tau / h;
for j = 1:n
    % Fictitious value from Neumann condition.
    Vfict = 2*h + V(m,j);
    Vj = [ V(:,j); Vfict ];
    % First row is zero by the Dirichlet condition.
    for i = 2:m+1 
        diff1 = (Vj(i+1) - Vj(i-1));
        diff2 = (Vj(i+1) - 2*Vj(i) + Vj(i-1));
        V(i,j+1) = Vj(i) ...
             + (lambda*sigma^2 * x(i)^2/2) * diff2  ...
             + (r*mu * x(i))/2 * diff1 - r*tau * Vj(i);
    end   
end
clf
for j = index_times
    str = sprintf("t = %.2f", t(j));
    plot(x, V(:, j), displayname=str) 
    hold on
end
title('Black-Scholes instability')
xlabel('stock price'),  ylabel('option value')
axis tight,  grid on
legend(location="northwest")
```

```{code-cell}
:tags: remove-output
clf
plot(x, V(:,1))
hold on,  grid on
axis([0, 8, 0, 6])
title('Black-Scholes solution...?') 
xlabel('stock price'),  ylabel('option value')
vid = VideoWriter("figures/black-scholes-8.mp4","MPEG-4");
vid.Quality = 85;
open(vid);
for frame = 1:10:n+1
    cla, plot(x, V(:, frame))
    str = sprintf("t = %.2f", t(frame));
    text(0.4, 5.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

![Trouble in Black–Scholes solution](figures/black-scholes-8.mp4)

This so-called solution is nonsense!
``````

### 11.2 @section-diffusion-methodlines

(demo-methodlines-heatFE-matlab)=
``````{dropdown} @demo-methodlines-heatFE
:open:
Let's implement the method of {numref}`Example {number} <example-methodlines-heatFE>` with second-order space semidiscretization.

```{code-cell}
m = 100;  
[x, Dx, Dxx] = diffper(m, [0, 1]);
Ix = eye(m);
```

Next we set an initial condition. It isn't mathematically periodic, but the end values and derivatives are so small that for numerical purposes it may as well be.

```{code-cell}
tfinal = 0.15;  n = 2400;  
tau = tfinal / n;  t = tau * (0:n)';
U = zeros(m, n+1);
U(:, 1) = exp( -60*(x - 0.5).^2 );
```

The Euler time stepping simply multiplies the solution vector by the constant matrix in {eq}`Eulerxx` at each time step. Since that matrix is sparse, we will declare it as such, even though the run-time savings may not be detectable for this small value of $m$.

```{code-cell}
A = sparse(Ix + tau * Dxx);
for j = 1:n
    U(:, j+1) = A * U(:,j);
end

index_times = 1:10:31;
show_times = t(index_times);
clf
for j = index_times
    str = sprintf("t = %.2e", t(j));
    plot(x, U(:, j), displayname=str) 
    hold on
end
legend(location="northwest")
xlabel('x'), ylabel('u(x,t)')
title('Heat equation by forward Euler')
```

You see above that things seem to start well, with the initial peak widening and shrinking. But then there is a nonphysical growth in the solution.

```{code-cell}
:tags: [remove-output, hide-input]
clf
index_times = 1:101;
plot(x, U(:, 1))
hold on,  grid on
axis([0, 1, -1, 2])
title('Heat equation by forward Euler') 
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/diffusionFE.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for frame = index_times
    cla, plot(x, U(:, frame))
    str = sprintf("t = %.3f", t(frame));
    text(0.05, 0.92, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```

![Instability in Euler solution](figures/diffusionFE.mp4)

The growth in norm is exponential in time.

```{code-cell}
M = max(abs(U), [], 1);     % max in each column
clf,  semilogy(t, M)
xlabel('t'), ylabel('max_x |u(x,t)|') 
title('Nonphysical growth')
```
``````

(demo-methodlines-heatBE-matlab)=
``````{dropdown} @demo-methodlines-heatBE
:open:
Now we apply backward Euler to the heat equation. Mathematically this means multiplying by the *inverse* of a matrix, but we interpret that numerically as a linear system solution. We will reuse the setup from @demo-methodlines-heatFE. 

```{code-cell}
B = sparse(Ix - tau * Dxx);
[l, u] = lu(B);
for j = 1:n
    U(:, j+1) = u \ (l \ U(:, j));
end

index_times = 1:600:n+1;
show_times = t(index_times);
clf
for j = index_times
    str = sprintf("t = %.2e", t(j));
    plot(x, U(:, j), displayname=str) 
    hold on
end
legend(location="northwest")
xlabel('x'), ylabel('u(x,t)')
title('Heat equation by backward Euler')
```

```{code-cell}
:tags: [hide-input, remove-output]
clf
index_times = 1:24:n+1;
plot(x, U(:, 1))
hold on,  grid on
axis([0, 1, -0.25, 1])
title('Heat equation by backward Euler') 
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/diffusionBE.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for frame = index_times
    cla, plot(x, U(:, frame))
    str = sprintf("t = %.3f", t(frame));
    text(0.05, 0.92, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid) 
```

![Stable Backward Euler solution](figures/diffusionBE.mp4)

This solution looks physically plausible, as the large concentration in the center diffuses outward until the solution is essentially constant. Observe that the solution remains periodic in space for all time.
``````
(demo-methodlines-auto-matlab)=
``````{dropdown} @demo-methodlines-auto
:open:
We set up the semidiscretization and initial condition in $x$ just as before.

```{code-cell}
m = 100;  
[x, Dx, Dxx] = diffper(m, [0, 1]);
Ix = eye(m);
u0 = exp( -60 * (x - 0.5).^2 );
```

Now, however, we apply a standard solver called `ode45` to the initial-value problem $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
tfinal = 0.05;
f = @(t, u) Dxx * u;
sol = ode45(f, [0, tfinal], u0);
u = @(t) deval(sol, t);

clf
for t = linspace(0, 0.05, 5)
    str = sprintf("t = %.3f", t);
    plot(x, u(t), displayname=str)
    hold on
end
xlabel("x"),  ylabel("u(x,t)")
legend()
title("Heat equation by ode45")
```

The solution appears to be correct. But the number of time steps that were selected automatically is surprisingly large, considering how smoothly the solution changes.

```{code-cell}
time_steps_ode45 = length(sol.x) - 1
```

Now we apply a different solver called `ode15s`.

```{code-cell}
sol = ode15s(f, [0, tfinal], u0);
u = @(t) deval(sol, t);
time_steps_ode15s = length(sol.x) - 1
```

The number of steps selected was reduced by a factor of 15!
``````

### 11.3 @section-diffusion-absstab
(demo-absstab-regions-matlab)=
``````{dropdown} @demo-absstab-regions
:open:
Euler and Backward Euler time-stepping methods were used to solve $\mathbf{u}'=\mathbf{D}_{xx}\mathbf{u}$.

```{code-cell}
m = 40;  
[x, Dx, Dxx] = diffper(m, [0, 1]);
```

The eigenvalues of this matrix are real and negative:

```{code-cell}
lambda = eig(Dxx);
clf
plot(real(lambda), imag(lambda), 'o')
axis equal, grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Eigenvalues')
```

The Euler method is absolutely stable in the region $|\zeta+1| \le 1$ in the complex plane:

```{code-cell}
:tags: hide-input
phi = linspace(0, 2*pi, 361);
z = exp(1i*phi) - 1;   % unit circle shifted to the left by 1
fill(real(z), imag(z), [.8, .8, 1])
axis equal,  grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Stability region')
```

In order to get inside this region, we have to find $\tau$ such that $\lambda \tau > -2$ for all eigenvalues $\lambda$. This is an upper bound on $\tau$. 

```{code-cell}
lambda_min = min(lambda)
max_tau = -2 / lambda_min
```

Here we plot the resulting values of $\zeta=\lambda \tau$. 

```{code-cell}
zeta = lambda * max_tau;
hold on
plot(real(zeta), imag(zeta), 'o')
title('Stability region and \zeta values')
```

In backward Euler, the region is $|\zeta-1|\ge 1$. Because they are all on the negative real axis, all of the $\zeta$ values will fit no matter what $\tau$ is chosen.

```{code-cell}
:tags: hide-input
clf,  fill([-6, 6, 6, -6],[-6, -6, 6, 6],[.8, .8, 1])
hold on
z = exp(1i*phi) + 1;   % unit circle shifted right by 1
fill(real(z), imag(z), 'w')
plot(real(zeta), imag(zeta), 'o')
axis equal,  axis([-4, 2, -3, 3]), grid on
xlabel('Re \lambda'),  ylabel('Im \lambda')
title('Stability region and \zeta values')
```
``````


### 11.4 @section-diffusion-stiffness
(demo-stiffness-oregon-matlab)=
``````{dropdown} @demo-stiffness-oregon
:open:
In {numref}`Example {number} <example-stiffness-oregon>` we derived a Jacobian matrix for the Oregonator model. Here is a numerical solution of the ODE.

```{code-cell}
:tags: hide-input
q = 8.375e-6;  s = 77.27;  w = 0.161;
f = @(t, u, p) [ s*(u(2)- u(1) * u(2) + u(1) - q * u(1)^2);...
    (-u(2) - u(1) * u(2) + u(3)) / s; ...
    w*(u(1) - u(3)) ];

sol = ode15s(f, [0, 500], [1; 2; 3]);
clf,  semilogy(sol.x, sol.y)
xlabel("t"),  ylabel("u(t)")
title("Oregonator solution")
```

At each value of the numerical solution, we can compute the eigenvalues of the Jacobian. Here we plot all of those eigenvalues in the complex plane.

```{code-cell}
:tags: hide-input
J = @(u) [ -s*(u(2)+1-2*q*u(1)), s*(1-u(1)), 0; ...
    -u(2)/s, (-1-u(1))/s, 1/s; ...
    w,0,-w];

t = sol.x;
u = sol.y;
lambda = zeros(length(t) - 1, 3);
for j = 1:length(t)-1
    lambda(j, :) = eig(J(u(:, j)));
end
plot3(real(lambda), imag(lambda), t(1:end-1), 'o')
xlabel("Re \lambda"),  ylabel("Im \lambda"),  zlabel("t")
title("Oregonator eigenvalues")
```

You can see that there is one eigenvalue that ranges over a wide portion of the negative real axis and dominates stability considerations.
``````

(demo-stiffness-explicit-matlab)=
``````{dropdown} @demo-stiffness-explicit
:open:
The `ode15s` solver is good for stiff problems and needs few time steps to solve the Oregonator from @demo-stiffness-oregon.

```{code-cell}
tic,  sol = ode15s(f, [0, 26], [1; 2; 3]); 
time_ode15s = toc
num_steps_ode15s = length(sol.x) - 1
```

But if we apply {numref}`Function {number} <function-rk23>` to the problem, the step size will be made small enough to cope with the large negative eigenvalue. 

```{code-cell}
ivp.ODEFcn = @(t, u, p) f(t, u);
ivp.InitialTime = 0;
ivp.InitialValue = [1; 2; 3];
ivp.Parameters = [];
tic, [t, u] = rk23(ivp, 0, 26, 1e-5);
time_rk23 = toc
num_steps_rk23 = length(t) - 1
```

Starting from the eigenvalues of the Jacobian matrix, we can find an effective $\zeta(t)$ by multiplying with the local time step size. The values of $\zeta(t)$ for each time level are plotted below and color coded by component of the diagonalized system.

```{code-cell}
:tags: hide-input
zeta = zeros(length(t) - 1, 3);
for j = 1:length(t)-1
    lambda = eig(J(u(:, j)));
    zeta(j, :) = (t(j+1) - t(j)) * lambda;
end
plot(zeta, 'o')
axis equal, grid on
xlabel('Re \zeta'),  ylabel('Im \zeta')
title("Oregonator stability")
```

Roughly speaking, the $\zeta$ values stay within or close to the RK2 stability region in {numref}`figure-stabreg_bd_rk`. Momentary departures from the region are possible, but time stepping repeatedly in that situation would cause instability. 
``````

### 11.5 @section-diffusion-boundaries
(demo-boundaries-heat-matlab)=
``````{dropdown} @demo-boundaries-heat
:open:
First, we define functions for the PDE and each boundary condition.

```{code-cell}
phi = @(t, x, u, ux, uxx) uxx;
ga = @(u, ux) u;
gb = @(u, ux) u - 2;
```

Our next step is to write a function to define the initial condition. This one satisfies the boundary conditions exactly.

```{code-cell}
init = @(x) 1 + sin(pi * x/2) + 3 * (1 - x.^2) .* exp(-4*x.^2);
```

Now we can use {numref}`Function {number} <function-parabolic>` to solve the problem.

```{code-cell}
[x, u] = parabolic(phi, [-1, 1], 60, ga, gb, [0, 0.75], init);

clf
for t = 0:0.1:0.5
    str = sprintf("t = %.2f", t);
    plot(x, u(t), displayname=str)
    hold on
end
xlabel("x"),  ylabel("u(x,t)")
legend()
title("Heat equation with Dirichlet boundaries")
```

```{code-cell}
:tags: remove-output
clf
plot(x, u(0))
hold on,  grid on
axis([-1, 1, 0, 4.2])
title('Heat equation with Dirichlet boundaries') 
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/boundaries-heat.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 0.75, 201)
    cla, plot(x, u(t))
    str = sprintf("t = %.3f", t);
    text(-0.9, 3.8, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid) 
```
![Heat equation with Dirichlet boundaries](figures/boundaries-heat.mp4)

``````


(demo-boundaries-bratu-matlab)=
``````{dropdown} @demo-boundaries-bratu
:open:

```{code-cell}
phi = @(t, x, u, ux, uxx) u.^2 + uxx;
ga = @(u, ux) u;
gb = @(u, ux) ux;

init = @(x) 400 * x.^4 .* (1 - x).^2;
[x, u] = parabolic(phi, [0, 1], 60, ga, gb, [0, 0.1], init);
```

```{code-cell}
:tags: remove-output
clf
plot(x, u(0))
hold on,  grid on
axis([0, 1, 0, 10])
title("Heat equation with source")
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/boundaries-source.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 0.1, 101)
    cla, plot(x, u(t))
    str = sprintf("t = %.3f", t);
    text(0.05, 9.2, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid) 
```

![Heat equation with source](figures/boundaries-source.mp4)
``````

(demo-boundaries-bs-matlab)=
``````{dropdown} @demo-boundaries-bs
:open:

```{code-cell}
K = 3;  sigma = 0.06;  r = 0.08;  Smax = 8;
phi = @(t, x, u, ux, uxx) sigma.^2/2 * (x.^2 .* uxx) + r * x.*ux - r*u;
ga = @(u, ux) u;
gb = @(u, ux) ux - 1;
```

```{code-cell}
init = @(x) max(0, x - K);
[x, u] = parabolic(phi, [0, Smax], 80, ga, gb, [0, 15], init);
```

```{code-cell}
:tags: remove-output
clf
plot(x, u(0))
hold on,  grid on
axis([0, Smax, -0.1, 8])
title("Black–Scholes equation with boundaries")
xlabel('x'),  ylabel('u(x,t)')
vid = VideoWriter("figures/boundaries-bs.mp4", "MPEG-4");
vid.Quality = 85;
open(vid);
for t = linspace(0, 15, 151)
    cla, plot(x, u(t))
    str = sprintf("t = %.1f", t);
    text(0.5, 7.1, str);
    writeVideo(vid, frame2im(getframe(gcf)));
end
close(vid)
```
![Black–Scholes equation with boundaries](figures/boundaries-bs.mp4)

Recall that $u$ is the value of the call option, and time runs backward from the strike time. The longer the horizon, the more value the option has due to anticipated growth in the stock price.
``````