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
[**Demo %s**](#demo-condition-roots)


The polynomial $p(x) = \frac{1}{3}(x-1)(x-1-\epsilon)$ has roots $1$ and $1+\epsilon$. For small values of $\epsilon$, the roots are ill-conditioned. 


```{code-cell} matlab
ep = 1e-6;
a = 1/3;             % coefficients of p...
b = (-2 - ep) / 3;   % ...
c = (1 + ep) / 3;    % ...in ascending order
```

Here are the roots as computed by the quadratic formula.

```{code-cell}
d = sqrt(b^2 - 4*a*c);
format long   % show all digits
r1 = (-b - d) / (2*a)
r2 = (-b + d) / (2*a)
```

The display of `r2` suggests that the last five digits or so are inaccurate. The relative errors are
```{tip}
:class: dropdown
Putting values inside square brackets creates a vector.
```

```{code-cell}
format short e
err = abs(r1 - 1) ./ abs(1)
err = abs(r2 - (1 + ep)) ./ abs(1 + ep)
```

The condition number of each root is 
$$
\kappa(r_i) = \frac{|r_i|}{|r_1-r_2|} \approx \frac{1}{\epsilon}. 
$$
Thus, relative error in the data at the level of roundoff can grow in the result to be roughly

```{code-cell}
eps / ep
```

This matches the observation pretty well.
