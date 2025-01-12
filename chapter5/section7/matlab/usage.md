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
[**Demo %s**](#demo-adapt-usage)

We'll integrate the function from {numref}`Demo %s <demo-adapt-motive>`.

```{code-cell}
f = @(x) (x + 1).^2 .* cos((2 * x + 1) ./ (x - 4.3));
```

We perform the integration and show the nodes selected underneath the curve.

```{code-cell}
[Q, t] = intadapt(f, 0, 4, 0.001);
clf, fplot(f, [0, 4], 2000)
hold on
stem(t, f(t), '.-')
title('Adaptive node selection')
xlabel('x'), ylabel('f(x)')
fprintf("number of nodes = %d", length(t))
```

The error turns out to be a bit more than we requested. It's only an estimate, not a guarantee.

```{code-cell}
I = integral(f, 0, 4, abstol=1e-14, reltol=1e-14);    % 'exact' value
fprintf("error = %.2e", abs(Q - I))
```

Let's see how the number of integrand evaluations and the error vary with the requested tolerance.

```{code-cell}
tol = 1 ./ 10.^(4:14)';
err = zeros(size(tol));
n = zeros(size(tol));
for i = 1:length(tol)
    [A, t] = intadapt(f, 0, 4, tol(i));
    err(i) =  I - A;
    n(i) = length(t);
end
disp(table(tol, err, n, variableNames=["tolerance", "error", "number of nodes"]))
```

As you can see, even though the errors are not smaller than the estimates, the two columns decrease in tandem. If we consider now the convergence not in $h$, which is poorly defined now, but in the number of nodes actually chosen, we come close to the fourth-order accuracy of the underlying Simpson scheme.

```{code-cell}
clf
loglog(n, abs(err), "-o", displayname="results")
xlabel("number of nodes"), ylabel("error")
title("Convergence of adaptive integration")
order4 = 0.1 * abs(err(end)) * (n / n(end)).^(-4);
hold on
loglog(n, order4, "k--", displayname="O(n^{-4})")
legend();
```
