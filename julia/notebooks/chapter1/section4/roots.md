---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---

For this example we will use the `Polynomials` package, which is installed by the `FNC` package.  

```{tip}
In the rest of the book, we do not show the `using` statement needed to load the book's package, but you will need to enter it if you want to run the codes yourself.
```

Our first step is to construct a polynomial with six known roots.

```{code-cell}
using Polynomials
r = [-2.0, -1, 1, 1, 3, 6]
p = fromroots(r)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
r̃ = sort(roots(p))   # type r\tilde and then press Tab
```

```{index} ! Julia; @., Julia; broadcasting
```

Here are the relative errors in each of the computed roots. 

```{tip}
The `@.` notation at the start means to do the given operations on each element of the given vectors.
```

```{code-cell}
println("Root errors:") 
@. abs(r - r̃) / r
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply `fromroots` to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
p̃ = fromroots(r̃)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell}
c,c̃ = coeffs(p), coeffs(p̃)
println("Coefficient errors:") 
@. abs(c - c̃) / c
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
