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
[**Demo %s**](#demo-float-double)


Python has native `int` and `float` types.

```{code-cell} ipython3
print(f"The type of {1} is {type(1)}")
print(f"The type of {float(1)} is {type(1.0)}")
```

The `numpy` package has its own `float` types:

```{code-cell} ipython3
one = float64(1)
print(f"The type of {one} is {type(one)}")
```

Both `float` and `float64` are double precision, using 64 binary bits per value. Although it is not normally necessary to do so, we can deconstruct a float into its significand and exponent:

```{code-cell} ipython3
x = 3.14
mantissa, exponent = frexp(x)
print(f"significand: {mantissa * 2}, exponent: {exponent - 1}")
```

```{code-cell} ipython3
mantissa, exponent = frexp(x / 8)
print(f"significand: {mantissa * 2}, exponent: {exponent - 1}")
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon, given here for double precision:

```{code-cell} ipython3
mach_eps = finfo(float).eps
print(f"machine epsilon is {mach_eps:.4e}")
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell} ipython3
print(f"machine epsilon is 2 to the power {log2(mach_eps)}")
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the significand. The actual range of positive values in double precision is

```{code-cell}
finf = finfo(float)
print(f"range of positive values: [{finf.tiny}, {finf.max}]")
```

For the most part you can mix integers and floating-point values and get what you expect.

```{code-cell}
1/7
```

```{code-cell}
37.3 + 1
```

```{code-cell}
2**(-4)
```

You can convert a floating value to an integer by wrapping it in `int`.

```{code-cell}
int(3.14)
```
