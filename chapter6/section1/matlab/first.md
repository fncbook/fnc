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
[**Demo %s**](#demo-basics-first)


Let's use it to define and solve an initial-value problem for $u'=\sin[(u+t)^2]$ over $t \in [0,4]$, such that $u(0)=-1$. To create an initial-value problem for $u(t)$, you must create an `ode` with a function that computes $u'$ and an initial condition for $u$. Then you create a solution by calling `solve` with a time interval. 
```{tip}
:class: dropdown
Most real ODE problems contain parameters that are constant during the solution but that can change from one problem instance to the next. Accordingly, we define the ODE function below to accept a third argument, `p`, which is a vector of parameters. We always include this argument for consistency, even when there are no parameters.
```

```{index} ! MATLAB; ode, ! MATLAB; solve
```

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u, p) sin((t + u)^2);
ivp.InitialTime = 0;
ivp.InitialValue = -1;
sol = solve(ivp, 0, 4);
```

The resulting solution object has fields `Time` and `Solution` that contain the approximate values of the solution at automatically chosen times in the interval you provided.

```{code-cell}
clf
plot(sol.Time, sol.Solution, '-o')
xlabel("t")
ylabel("u(t)")
title(("Solution of an IVP"));
```

You might want to know the solution at particular times other than the ones selected by the solver. That requires an interpolation, which is done by `solutionFcn`.

```{code-cell}
u = solutionFcn(ivp, 0, 10);
u(0:5)
```
