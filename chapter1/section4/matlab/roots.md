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
[**Demo %s**](#demo-stability-roots)



Our first step is to construct a polynomial with six known roots.
```{tip}
:class: dropdown
The `'` operator is used for transposition. Here, we want to make `r` a column vector.
```

```{code-cell}
r = [-2 ,-1, 1, 1, 3, 6]';
p = poly(r)
```

Now we use a standard numerical method for finding those roots, pretending that we don't know them already. This corresponds to $\tilde{y}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
rr = sort(roots(p))   
```

Here are the relative errors in each of the computed roots. 
```{tip}
:class: dropdown
The `./` operator is used for element-wise division.
```

```{code-cell}
disp("Root errors:") 
abs(r - rr) ./ r
```

It seems that the forward error is acceptably close to machine epsilon for double precision in all cases except the double root at $x=1$. This is not a surprise, though, given the poor conditioning at such roots.

Let's consider the backward error. The data in the rootfinding problem is the polynomial coefficients. We can apply `poly` to find the coefficients of the polynomial (that is, the data) whose roots were actually computed by the numerical algorithm. This corresponds to $\tilde{x}$ in {numref}`Definition {number} <definition-stability-backward>`. 

```{code-cell}
pp = poly(rr)
```

We find that in a relative sense, these coefficients are very close to those of the original, exact polynomial:

```{code-cell}
disp("Coefficient errors:") 
abs(p - pp) ./ abs(p)
```

In summary, even though there are some computed roots relatively far from their correct values, they are nevertheless the roots of a polynomial that is very close to the original.
