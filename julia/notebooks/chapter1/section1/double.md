---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
In Julia, `1` and `1.0` are different values, because they have different types:

```{code-cell}
@show typeof(1);
@show typeof(1.0);
```

The standard choice for floating-point values is `Float64`, which is double precision using 64 binary bits. We can see all the bits by using `bitstring`.

```{code-cell}
bitstring(1.0)
```


The first bit determines the sign of the number:

```{tip}
Square brackets concatenate the contained values into vectors.
```

```{code-cell}
[bitstring(1.0), bitstring(-1.0)]
```

The next 11 bits determine the exponent (scaling) of the number, and so on.

```{code-cell}
[bitstring(1.0), bitstring(2.0)]
```

The sign bit, exponent, and significand in {eq}`floatpoint` are all directly accessible.

```{code-cell}
x = 3.14
@show sign(x), exponent(x), significand(x);
```

```{code-cell}
x = x / 8
@show sign(x), exponent(x), significand(x);
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon. You can get its value from the `eps` function in Julia. By default, it returns the value for double precision.

```{tip}
To call a function, including `eps`, you must use parentheses notation, even when there are no input arguments.
```

```{code-cell} julia
eps()
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell}
log2(eps())
```

The spacing between adjacent floating-point values is proportional to the magnitude of the value itself. This is how relative precision is kept roughly constant throughout the range of values. You can get the adjusted spacing by calling `eps` with a value.

```{code-cell}
eps(1.618)
```

```{code-cell}
eps(161.8)
```

```{code-cell}
nextfloat(161.8)
```

```{code-cell}
ans - 161.8
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the mantissa. The actual range of positive values in double precision is

```{code-cell}
@show floatmin(), floatmax();
```

For the most part you can mix integers and floating-point values and get what you expect.

```{code-cell}
1/7
```

```{code-cell}
37.3 + 1
```

```{code-cell}
2^(-4)
```

There are some exceptions. A floating-point value can't be used as an index into an array, for example, even if it is numerically equal to an integer. In such cases you use `Int` to convert it.

```{code-cell}
@show 5.0, Int(5.0);
```

If you try to convert a noninteger floating-point value into an integer you get an `InexactValue` error. This occurs whenever you try to force a type conversion that doesn't make clear sense.
