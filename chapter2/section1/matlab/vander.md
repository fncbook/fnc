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
[**Demo %s**](#demo-interp-vander)

We create two column vectors for data about the population of China. The first has the years of census data and the other has the population, in millions of people.

```{code-cell}
year = [1982; 2000; 2010; 2015]; 
pop = [1008.18; 1262.64; 1337.82; 1374.62];
```

It's convenient to measure time in years since 1980. 

```{code-cell}
t = year - 1980;
y = pop;
```

```{index} ! MATLAB; vander
```

Now we have four data points $(t_1,y_1),\dots,(t_4,y_4)$, so $n=4$ and we seek an interpolating cubic polynomial. We construct the associated Vandermonde matrix:

```{code-cell}
V = vander(t)
```

:::{index} ! MATLAB; \\
:::

To solve for the vector of polynomial coefficients, we use a backslash to solve the linear system:
```{tip}
:class: dropdown
A **backslash** `\` is used to solve a linear system of equations.
```

```{code-cell}
c = V \ y
```

The algorithms used by the backslash operator are the main topic of this chapter. As a check on the solution, we can compute the *residual*.

```{code-cell}
y - V * c
```

Using floating-point arithmetic, it is not realistic to expect exact equality of quantities; a relative difference comparable to $\macheps$ is all we can look for.

By our definitions, the elements of `c` are coefficients in descending-degree order for the interpolating polynomial. We can use the polynomial to estimate the population of China in 2005:

```{code-cell}
p = @(t) polyval(c, t - 1980);  % include the 1980 time shift
p(2005)
```

The official population value for 2005 was 1303.72, so our result is rather good. 

:::{index} ! MATLAB; scatter
:::

We can visualize the interpolation process. First, we plot the data as points.
```{tip}
:class: dropdown
The `scatter` function creates a scatter plot of points; you can specify a line connecting the points as well.
```

```{code-cell}
scatter(year, y)
xlabel("years since 1980")
ylabel("population (millions)")
title(("Population of China"));
```

:::{index} ! MATLAB; linspace
:::

We want to superimpose a plot of the polynomial. We do that by evaluating it at a vector of points in the interval.
```{tip}
:class: dropdown
The `linspace` function constructs evenly spaced values given the endpoints and the number of values.
```

```{code-cell}
tt = linspace(1980, 2015, 500);    % 500 times in the interval [1980, 2015]
yy = p(tt);                        % evaluate p at all the vector elements
yy(1:4)
```

```{index} ! MATLAB; hold on, ! MATLAB; plot
```

Now we use `plot!` to add to the current plot, rather than replacing it.
```{tip}
:class: dropdown
Use `hold on` to add to an existing plot rather than replacing it.
The `plot` function plots lines connecting the given $x$ and $y$ values; you can also specify markers at the points.
```

```{code-cell}
hold on 
plot(tt, yy)
legend("data", "interpolant", "location", "northwest");
```
