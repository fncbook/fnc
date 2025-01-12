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
[**Demo %s**](#demo-systems-predator)

We encode the predator–prey equations via a function, defined here externally.

```{literalinclude} f63_predprey.m
:language: matlab
```

The values of `alpha` and `beta` are parameters that influence the solution of the IVP. We use the `Parameters` field of the IVP object to define them for the solver, which in turn passes them as the third argument into our ODE function. 

```{code-cell}
u0 = [1; 0.01];    % column vector
p = [0.1, 0.25];
ivp = ode;
ivp.ODEFcn = @f63_predprey;
ivp.InitialValue = u0;
ivp.Parameters = p;
sol = solve(ivp, 0, 60);
size(sol.Solution)
```

Each column of the `Solution` field is the solution vector $\mathbf{u}$ at a particular time; each row is a component of $\mathbf{u}$ over all time.

```{code-cell}
clf
plot(sol.Time, sol.Solution)
xlabel("t")
ylabel("u(t)")
title('Predator-prey solution')
legend('prey', 'predator');
```

We can also use {numref}`Function {number} <function-euler>` to find the solution.

```{code-cell}
[t, u] = eulerivp(ivp, 0, 60, 1200);
```

```{code-cell}
hold on
plot(t, u, '.')
```

Notice above that the accuracy of the Euler solution deteriorates rapidly.

When there are just two components, it's common to plot the solution in the _phase plane_, i.e., with $u_1$ and $u_2$ along the axes and time as a parameterization of the curve.

```{code-cell}
clf
plot(u(1, :), u(2, :))
title("Predator-prey in the phase plane")
xlabel("y")
ylabel(("z"));
```

From this plot we can deduce that the solution approaches a periodic one, which in the phase plane is represented by a closed loop.
