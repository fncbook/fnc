---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
```{code-cell}
:tags: [remove-cell]
exec(open("../../../python/FNC_init.py").read())
```
[**Demo %s**](#demo-interp-vander)

We create two vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell} 
year = arange(1980, 2020, 10)   # from 1980 to 2020 by 10
pop = array([984.736, 1148.364, 1263.638, 1330.141])
```

It's convenient to measure time in years since 1980. 

```{code-cell} 
t = year - 1980
y = pop
```

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix: 

```{code-cell} 
V = vander(t)
print(V)
```

To solve a linear system $\mathbf{V} \mathbf{c} = \mathbf{y}$ for the vector of polynomial coefficients, we use `solve` (imported from `numpy.linalg`):

```{code-cell} 
c = linalg.solve(V, y)
print(c)
```
The algorithms used by `solve` are the main topic of this chapter. As a check on the solution, we can compute the *residual* $\mathbf{y} - \mathbf{V} \mathbf{c}$, which should be small (near machine precision).

```{tip}
:class: dropdown
Matrix multiplication in NumPy is done with `@` or `matmul`.
```

```{code-cell} 
print(y - V @ c)
```

By our definitions, the coefficients in `c` are given in descending order of power in $t$. We can use the resulting polynomial to estimate the population of China in 2005:

```{code-cell} 
p = poly1d(c)          # construct a polynomial
print(p(2005 - 1980))     # apply the 1980 time shift
```

The official figure was 1303.72, so our result is rather good.

We can visualize the interpolation process. First, we plot the data as points. Then we add a plot of the interpolant, taking care to shift the $t$ variable back to actual years.

```{code-cell} 
scatter(year, y, color="k", label="data");
tt = linspace(0, 30, 300)   # 300 times from 1980 to 2010
plot(1980 + tt, p(tt), label="interpolant");
xlabel("year");
ylabel("population (millions)");
title("Population of China");
legend();
```
