---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

```{code-cell}
using Plots, SpecialFunctions
J₃(x) = besselj(3, x)
plot(J₃, 0, 20;
    title="Bessel function",
    xaxis=(L"x"),  yaxis=(L"J_3(x)"),  grid=:xy)
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `nlsolve` from the `NLsolve` package to find these roots accurately. It uses vector variables, so we have to code accordingly.
```{tip}
Type `\omega` followed by <kbd>Tab</kbd> to get the character `ω`.
The argument `ftol=1e-14` below is called a **keyword argument**. Here it sets a goal for the maximum value of $|f(x)|$.
```

```{code-cell}
using NLsolve
ω = []
for guess = [6., 10. ,13., 16., 19.]
    s = nlsolve(x -> J₃(x[1]), [guess], ftol=1e-14)
    append!(ω, s.zero)
end
```

```{code-cell}
y = J₃.(ω)
@pt :header=["root estimate", "function value"] [ω y]
```

```{code-cell}
scatter!(ω, y, title="Bessel function with roots")
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
r = []
for guess = [3., 6., 10., 13.]
    f(x) = J₃(x[1]) - 0.2
    s = nlsolve(f, [guess], ftol=1e-14)
    append!(r, s.zero)
end
scatter!(r, J₃.(r), title="Roots and other Bessel values")
```
