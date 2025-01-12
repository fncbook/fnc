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
[**Demo %s**](#demo-fitting-pirate)

```{code-cell}
k = (1:100)';
a = 1./k.^2;      % sequence
s = cumsum(a);    % cumulative summation
p = sqrt(6*s);
clf
plot(k, p, 'o-')
xlabel('k'), ylabel('p_k')
title(('Sequence converging to \pi'));
```

This graph suggests that maybe $p_k\to \pi$, but it's far from clear how close the sequence gets. It's more informative to plot the sequence of errors, $\epsilon_k= |\pi-p_k|$. By plotting the error sequence on a log-log scale, we can see a nearly linear relationship.

```{code-cell}
ep = abs(pi - p);    % error sequence
loglog(k, ep, 'o')
title('Convergence')
xlabel('k'), ylabel('|p_k - \pi|'), axis tight    
```

The straight line on the log-log scale suggests a power-law relationship where $\epsilon_k\approx a k^b$, or $\log \epsilon_k \approx b (\log k) + \log a$.

```{code-cell}
V = [ k.^0, log(k) ];    % fitting matrix
c = V \ log(ep)          % coefficients of linear fit
```

In terms of the parameters $a$ and $b$ used above, we have

```{code-cell}
a = exp(c(1)),  b = c(2)
```

It's tempting to conjecture that the slope $b\to -1$ asymptotically. Here is how the numerical fit compares to the original convergence curve.

```{code-cell}
hold on
loglog(k, a * k.^b)
legend('sequence', 'power-law fit');
```
