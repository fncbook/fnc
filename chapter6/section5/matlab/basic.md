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
[**Demo %s**](#demo-adapt-basic)

Let's run adaptive RK on  $u'=e^{t-u\sin u}$.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u, p) exp(t - u * sin(u));
ivp.InitialValue = 0;
a = 0;  b = 5;

[t, u] = rk23(ivp, a, b, 1e-5);
clf, plot(t, u)
xlabel("t");  ylabel("u(t)")
title(("Adaptive IVP solution"));
```

The solution makes a very abrupt change near $t=2.4$. The resulting time steps vary over three orders of magnitude.

```{code-cell}
Delta_t = diff(t);
semilogy(t(1:end-1), Delta_t) 
xlabel("t");  ylabel("step size")
title(("Adaptive step sizes"));
```

If we had to run with a uniform step size to get this accuracy, it would be

```{code-cell}
fprintf("minimum step size = %.2e", min(Delta_t))
```

On the other hand, the average step size that was actually taken was

```{code-cell}
fprintf("average step size = %.2e", mean(Delta_t))
```

We took fewer steps by a factor of almost 1000! Even accounting for the extra stage per step and the occasional rejected step, the savings are clear.

