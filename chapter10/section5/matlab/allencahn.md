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
[**Demo %s**](#demo-nonlinear-allencahn)


```{code-cell}
epsilon = 0.05;
phi = @(x, u, du_dx) (u.^3 - u) / epsilon;
ga = @(u, du) du;
gb = @(u, du) u - 1;
```

Finding a solution is easy at larger values of $\epsilon$.

```{code-cell}
init = linspace(-1, 1, 141)';
[x, u1] = bvp(phi, 0, 1, ga, gb, init);
clf,  plot(x, u1, displayname="\epsilon = 0.05")
xlabel('x'),  ylabel('u(x)')
title('Allen-Cahn solution') 
legend(location="northwest") 
```

However, finding a good initialization is not trivial for smaller values of $\epsilon$. Note below that the iteration stops without converging to a solution.

```{code-cell}
epsilon = 0.002;
phi = @(x, u, du_dx) (u.^3 - u) / epsilon;
[x, z] = bvp(phi, 0, 1, ga, gb, init);
```

The iteration succeeds if we use the first solution instead as the initialization here.

```{code-cell}
[x, u2] = bvp(phi, 0, 1, ga, gb, u1);
hold on,  plot(x, u2, displayname="\epsilon = 0.002")
```

In this case we can continue further.

```{code-cell}
epsilon = 0.0005;
phi = @(x, u, du_dx) (u.^3 - u) / epsilon;
[x, u3] = bvp(phi, 0, 1, ga, gb, u2);
plot(x, u3, displayname="\epsilon = 0.0005")

```
