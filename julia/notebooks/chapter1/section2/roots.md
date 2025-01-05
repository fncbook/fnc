---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 
```{tip}
The statement `x,y = 10,20` makes individual assignments to both `x` and `y`.
```


```{code-cell}
ϵ = 1e-6   # type \epsilon and then press Tab
a,b,c = 1/3,(-2-ϵ)/3,(1+ϵ)/3   # coefficients of p
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b^2 - 4a*c)
r₁ = (-b - d) / (2a)   # type r\_1 and then press Tab
r₂ = (-b + d) / (2a)
(r₁, r₂)
```

The relative errors in these values are 

```{code-cell}
@show abs(r₁ - 1) / abs(1);
@show abs(r₂ - (1+ϵ)) / abs(1+ϵ);
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly


```{code-cell}
eps() / ϵ
```

This matches the observation pretty well.
