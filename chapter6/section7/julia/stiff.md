---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-implicit-stiff)

The following simple ODE uncovers a surprise.

```{code-cell}
ivp = ODEProblem((u, p, t) -> u^2 - u^3, 0.005, (0, 400.0))
```

We will solve the problem first with the implicit AM2 method using $n=200$ steps.

```{code-cell}
tI, uI = FNC.am2(ivp, 200)

plot(tI, uI;
    label="AM2", legend=:bottomright,
    xlabel=L"t",  ylabel=L"u(t)")
```

Now we repeat the process using the explicit AB4 method.

```{code-cell}
tE, uE = FNC.ab4(ivp, 200)
scatter!(tE, uE, m=3, label="AB4", ylim=[-4, 2])
```

Once the solution starts to take off, the AB4 result goes catastrophically wrong.

```{code-cell}
uE[105:111]
```

We hope that AB4 will converge in the limit $h\to 0$, so let's try using more steps.

```{code-cell}
plt = scatter(tI, uI;
    m=3,  label="AM2, n=200",  legend=:bottomright,
    xlabel=L"t",  ylabel=L"u(t)")

for n in [1000, 1600]
    tE, uE = FNC.ab4(ivp, n)
    plot!(tE, uE, label="AM4, n=$n")
end
plt
```

So AB4, which is supposed to be _more_ accurate than AM2, actually needs something like 8 times as many steps to get a reasonable-looking answer!
