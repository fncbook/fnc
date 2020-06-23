---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Julia (fast start)
  language: julia
  name: julia-fast
---

# Polynomial roots

For this example we will use a publicly available package for working with polynomials. It should be available using the following line, if you have followed installation instructions for these scripts.

```{code-cell}
using Polynomials
```

Our first step is to construct a polynomial with six known roots.

```{code-cell}
r = [-2.0,-1,1,1,3,6]
@show p = fromroots(r);
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already.

```{code-cell}
@show r_computed = sort(roots(p));
```

Here are the relative errors in each of the computed roots. The `@.` notation at the start means essentially to do the given operations on each element of the given vectors.

```{code-cell}
@. abs(r - r_computed) / r
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem are the polynomial coefficients. We can apply `fromroots` to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm.

```{code-cell}
@show p_computed = fromroots(r_computed);
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell}
cp = coeffs(p)
cpc = coeffs(p_computed)
@. abs(cp-cpc)/cp
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
