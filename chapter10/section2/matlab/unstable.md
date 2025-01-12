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
[**Demo %s**](#demo-shooting-unstable)


```{code-cell}
:tags: [raises-exception]
ga = @(u, du) u + 1;
gb = @(u, du) u;
clf
warning off
for lambda = 16:4:28
    phi = @(x, u, du_dx) lambda^2 * (u + 1);
    [x, u, du_dx] = shoot(phi, 0.0, 1.0, ga, gb, [-1; 0], 1e-8);
    plot(x, u, displayname=sprintf("lambda=%d", lambda))
    hold on
    xlabel('x'),  ylabel('u(x)')
    title('Shooting instability')
    legend(location="northwest");
end
```

The numerical solution fails at the largest value of $\lambda$ because the initial condition became infinite.
