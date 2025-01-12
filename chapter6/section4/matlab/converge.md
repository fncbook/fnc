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
[**Demo %s**](#demo-rk-converge)

We solve the IVP $u'=\sin[(u+t)^2]$ over $0\le t \le 4$, with $u(0)=-1$.

```{code-cell}
ivp = ode;
ivp.ODEFcn = @(t, u, p) sin((t + u)^2);
ivp.InitialValue = -1;
a = 0;  b = 4;
```

We use a built-in solver to construct an accurate approximation to the exact solution.

```{code-cell}
ivp.AbsoluteTolerance = 1e-13;
ivp.RelativeTolerance = 1e-13;
u_ref = solutionFcn(ivp, a, b);
```

Now we perform a convergence study of our two Runge–Kutta implementations.

```{code-cell}
n = round(2 * 10.^(0:0.5:3)');
err = zeros(length(n), 2);
for i = 1:length(n)
    [t, u] = ie2(ivp, a, b, n(i));
    err(i, 1) = norm(u_ref(t) - u, Inf);
    [t, u] = rk4(ivp, a, b, n(i));
    err(i, 2) = norm(u_ref(t) - u, Inf);
end

disp(table(n, err(:, 1), err(:, 2), variableNames=["n", "IE2 error", "RK4 error"]))
```

The amount of computational work at each time step is assumed to be proportional to the number of stages. Let's compare on an apples-to-apples basis by using the number of $f$-evaluations on the horizontal axis.

```{code-cell}
clf, loglog([2*n 4*n], err, '-o')
hold on
loglog(2*n, 1e-5 * (n / n(end)) .^ (-2), '--')
loglog(4*n, 1e-10 * (n / n(end)) .^ (-4), '--')
xlabel("f-evaluations");  ylabel("inf-norm error")
title("Convergence of RK methods")
legend("IE2", "RK4", "O(n^{-2})", "O(n^{-4})", "location", "southwest");
```

The fourth-order variant is more efficient in this problem over a wide range of accuracy.
