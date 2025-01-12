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
[**Demo %s**](#demo-fdconverge-order12)

Let's observe the convergence of the formulas in {numref}`Example {number} <example-fd-converge-FD11>` and {numref}`Example {number} <example-fd-converge-FD12>`, applied to the function $\sin(e^{x+1})$ at $x=0$.

```{code-cell}
f = @(x) sin(exp(x + 1));
exact_value = exp(1) * cos(exp(1))
```

We'll compute the formulas in parallel for a sequence of $h$ values.

```{code-cell}
h = 5 ./ 10.^(1:6)';
FD1 = zeros(size(h));
FD2 = zeros(size(h));
for i = 1:length(h)
    h_i = h(i);
    FD1(i) = (f(h_i) - f(0)    ) / h_i;
    FD2(i) = (f(h_i) - f(-h_i)) / (2*h_i);
end
disp(table(h, FD1, FD2))
```

All that's easy to see from this table is that FD2 appears to converge to the same result as FD1, but more rapidly. A table of errors is more informative.

```{code-cell}
err1 = abs(exact_value - FD1);
err2 = abs(exact_value - FD2);
disp(table(h, err1, err2, variableNames=["h", "error in FD1", "error in FD2"]))
```

In each row, $h$ is decreased by a factor of 10, so that the error is reduced by a factor of 10 in the first-order method and 100 in the second-order method.

A graphical comparison can be useful as well. On a log-log scale, the error should (as $h\to 0$) be a straight line whose slope is the order of accuracy. However, it's conventional in convergence plots to show $h$ _decreasing_ from left to right, which negates the slopes.

```{code-cell}
clf
loglog(h, abs([err1 err2]), "o-")
set(gca, "xdir", "reverse")
order1 = 0.1 * err1(end) * (h / h(end)) .^ 1;
order2 = 0.1 * err2(end) * (h / h(end)) .^ 2;
hold on
loglog(h, order1, "--", h, order2, "--")
xlabel("h");  ylabel("error")
title("Convergence of finite differences")
legend("FD1", "FD2", "O(h)", "O(h^2)");
```
