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
[**Demo %s**](#demo-fitting-tempinterp)

Here are 5-year averages of the worldwide temperature anomaly as compared to the 1951–1980 average (source: NASA).

```{code-cell}
year = arange(1955,2005,5)
y = array([ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
    0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ])

fig, ax = subplots()
ax.scatter(year, y, color="k", label="data")
xlabel("year")
ylabel("anomaly (degrees C)")
title("World temperature anomaly");
```

A polynomial interpolant can be used to fit the data. Here we build one using a Vandermonde matrix. First, though, we express time as decades since 1950, as it improves the condition number of the matrix.

```{code-cell}
t = (year - 1950) / 10
V = vander(t)
c = linalg.solve(V, y)
print(c)
```

```{index} Python; plotting functions
```

The coefficients in vector `c` are used to create a polynomial. Then we create a function that evaluates the polynomial after changing the time variable as we did for the Vandermonde matrix.

```{code-cell}
p = poly1d(c)    # convert to a polynomial
tt = linspace(1955, 2000, 500)
ax.plot(tt, p((tt - 1950) / 10), label="interpolant")
ax.legend();
fig
```

As you can see, the interpolant does represent the data, in a sense. However it's a crazy-looking curve for the application. Trying too hard to reproduce all the data exactly is known as _overfitting_.

