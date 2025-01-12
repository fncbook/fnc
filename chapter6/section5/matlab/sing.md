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
[**Demo %s**](#demo-adapt-sing)

In {numref}`Demo %s <demo-basics-sing>` we saw an IVP that appears to blow up in a finite amount of time. Because the solution increases so rapidly as it approaches the blowup, adaptive stepping is required even to get close.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u, p) (t + u)^2;
ivp.InitialValue = 1;
[t, u] = rk23(ivp, 0, 1, 1e-5);
```

In fact, the failure of the adaptivity gives a decent idea of when the singularity occurs.

```{code-cell}
clf, semilogy(t, u)
xlabel("t");  ylabel("u(t)")
title("Adaptive solution near a singularity")

tf = t(end);
xline(tf, "linestyle", "--")
text(tf, 1e5, sprintf(" t = %.6f ", tf))
```
