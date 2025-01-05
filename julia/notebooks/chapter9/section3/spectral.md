---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell} 
:tags: [remove-cell]
using Logging
disable_logging(Logging.Warn);
```
On the left, we use a log-log scale, which makes second-order algebraic convergence $O(n^{-4})$ a straight line. On the right, we use a log-linear scale, which makes spectral convergence $O(K^{-n})$ linear.

```{code-cell} 
:tags: [hide-input]
n = 20:20:400
algebraic = @. 100 / n^4
spectral = @. 10 * 0.85^n
plot(n, [algebraic spectral], layout=(1, 2), subplot=1,
    xaxis=(L"n", :log10),  yaxis=(:log10, (1e-15, 1)),
    label=["algebraic" "spectral"],  title="Log-log")
plot!(n, [algebraic spectral], subplot=2,
    xaxis=L"n",  yaxis=(:log10, (1e-15, 1)),
    label=["algebraic" "spectral"],  title="log-linear")
```
