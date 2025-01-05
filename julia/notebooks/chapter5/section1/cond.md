---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering: 
  headings: false
---
In {numref}`Demo %s <demo-interpolation-global>` and {numref}`Demo %s <demo-interpolation-pwise>` we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
n = 18
t = range(-1, 1, n+1)
y = [zeros(9); 1; zeros(n - 9)];  # data for 10th cardinal function

scatter(t, y, label="data")
```

```{code-cell}
ϕ = Spline1D(t, y)
plot!(x -> ϕ(x), -1, 1;
    label="spline",
    xlabel=L"x",  ylabel=L"\phi(x)",
    title="Piecewise cubic cardinal function")
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
scatter(t, y, label="data")

ϕ = Polynomials.fit(t, y, n)
plot!(x -> ϕ(x), -1, 1;
    label="polynomial",  legend=:top,
    xlabel=L"x",  ylabel=L"\phi(x)", 
    title="Polynomial cardinal function")
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
