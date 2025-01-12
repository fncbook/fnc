---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```
[**Demo %s**](#demo-stiffness-explicit)

The `ode15s` solver is good for stiff problems and needs few time steps to solve the Oregonator from {numref}`Demo {number} <demo-stiffness-oregon>`.

```{code-cell}
tic,  sol = solve(ivp, 0, 26); 
time_ode15s = toc
num_steps_ode15s = length(sol.Time) - 1
```

But if we apply {numref}`Function {number} <function-rk23>` to the problem, the step size will be made small enough to cope with the large negative eigenvalue. 

```{code-cell}
tic, [t, u] = rk23(ivp, 0, 26, 1e-5);
time_rk23 = toc
num_steps_rk23 = length(t) - 1
```

Starting from the eigenvalues of the Jacobian matrix, we can find an effective $\zeta(t)$ by multiplying with the local time step size. The values of $\zeta(t)$ for each time level are plotted below and color coded by component of the diagonalized system.

```{code-cell}
:tags: [hide-input]
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
