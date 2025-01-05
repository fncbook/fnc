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
[**Demo %s**](#demo-finitediffs-fd-weights)

We will estimate the derivative of $\cos(x^2)$ at $x=0.5$ using five nodes.

```{code-cell}
t = [0.35, 0.5, 0.57, 0.6, 0.75]   # nodes
f = x -> cos(x^2)
df_dx = x -> -2 * x * sin(x^2)
exact_value = df_dx(0.5)
```

We have to shift the nodes so that the point of estimation for the derivative is at $x=0$. (To subtract a scalar from a vector, we must use the `.-` operator.)

```{code-cell}
w = FNC.fdweights(t .- 0.5, 1)
```

The finite-difference formula is a dot product (i.e., inner product) between the vector of weights and the vector of function values at the nodes.

```{code-cell}
fd_value = dot(w, f.(t))
```

We can reproduce the weights in the finite-difference tables by using equally spaced nodes with $h=1$. For example, here is a one-sided formula at four nodes.

```{code-cell}
FNC.fdweights(0:3, 1)
```

```{index} ! Julia; Rational
```

By giving nodes of type `Rational`, we can get exact values instead.

```{code-cell}
FNC.fdweights(Rational.(0:3), 1)
```
