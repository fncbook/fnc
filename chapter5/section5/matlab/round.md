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
[**Demo %s**](#demo-fdconverge-round)

Let $f(x)=e^{-1.3x}$. We apply finite-difference formulas of first, second, and fourth order to estimate $f'(0)=-1.3$.

```{code-cell}
f = @(x) exp(-1.3 * x);
exact = -1.3;

h = 10 .^ (-(1:12))';
FD = zeros(length(h), 3);
for i = 1:length(h)
    h_i = h(i);
    nodes = h_i * (-2:2);
    vals = f(nodes);
    FD(i, 1) = dot([0      0 -1   1    0] / h_i, vals);
    FD(i, 2) = dot([0    -1/2 0 1/2    0] / h_i, vals);
    FD(i, 3) = dot([1/12 -2/3 0 2/3 -1/12] / h_i, vals);
end
format long
disp(table(h, FD(:, 1), FD(:, 2), FD(:, 3), variableNames=["h", "FD1", "FD2", "FD4"]))
```

They all seem to be converging to $-1.3$. The convergence plot reveals some interesting structure to the errors, though.

```{code-cell}
err = abs(FD - exact);
clf
loglog(h, err, "o-")
set(gca, "xdir", "reverse")
order1 = 0.1 * err(end, 1) * (h / h(end)) .^ (-1);
hold on
loglog(h, order1, "k--")
xlabel("h");  ylabel("error")
title("FD error with roundoff")
legend("FD1", "FD2", "FD4", "O(1/h)", "location", "northeast");
```

Again the graph is made so that $h$ decreases from left to right. The errors are dominated at first by truncation error, which decreases most rapidly for the fourth-order formula. However, increasing roundoff error eventually equals and then dominates the truncation error as $h$ continues to decrease. As the order of accuracy increases, the crossover point moves to the left (greater efficiency) and down (greater accuracy).
