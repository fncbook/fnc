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
[**Demo %s**](#demo-basics-sing)


The equation $u'=(u+t)^2$ gives us some trouble.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u, p) (t + u)^2;
ivp.InitialTime = 0;
ivp.InitialValue = 1;
sol = solve(ivp, 0, 1);
```

The warning message we received can mean that there is a bug in the formulation of the problem. But if everything has been done correctly, it suggests that the solution may not exist past the indicated time. This is a possibility in nonlinear ODEs.

```{code-cell}
clf
semilogy(sol.Time, sol.Solution)
xlabel("t")
ylabel("u(t)")
title(("Finite-time blowup"));
```
