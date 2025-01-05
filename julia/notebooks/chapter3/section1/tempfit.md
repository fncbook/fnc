---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Here are the 5-year temperature averages again.

```{code-cell}
year = 1955:5:2000
t = @. (year - 1950) / 10
temp = [ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
          0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ]
```

```{index} Julia; \\
```

The standard best-fit line results from using a linear polynomial that meets the least-squares criterion.
```{tip}
Backslash solves overdetermined linear systems in a least-squares sense.
```

```{code-cell}
V = [ t.^0 t ]    # Vandermonde-ish matrix
@show size(V)
c = V \ temp
p = Polynomial(c)
```

```{code-cell}
f = yr -> p((yr - 1955) / 10)
scatter(year, temp, label="data",
    xlabel="year", ylabel="anomaly (degrees C)", leg=:bottomright)
plot!(f, 1955, 2000, label="linear fit")
```

If we use a global cubic polynomial, the points are fit more closely.

```{code-cell}
V = [ t[i]^j for i in 1:length(t), j in 0:3 ]   
@show size(V);
```

Now we solve the new least-squares problem to redefine the fitting polynomial.
```{tip}
The definition of `f` above is in terms of `p`. When `p` is changed, then `f` calls the new version.
```

```{code-cell}
p = Polynomial( V \ temp )
plot!(f, 1955, 2000, label="cubic fit")
```

If we were to continue increasing the degree of the polynomial, the residual at the data points would get smaller, but overfitting would increase.
