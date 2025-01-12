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
[**Demo %s**](#demo-float-double)

In MATLAB, values are double-precision floats unless declared otherwise. 

```{code-cell}
fprintf('1 has type: %s', class(1))
fprintf('1.0 has type: %s', class(1.0))
```

The spacing between floating-point values in $[2^n,2^{n+1})$ is $2^n \epsilon_\text{mach}$, where $\epsilon_\text{mach}$ is machine epsilon. Its value is predefined as `eps`.

```{tip}
:class: dropdown
While you can assign a different value to `eps`, doing so does not change any arithmetic. It's generally a bad idea. 
```

```{code-cell} 
eps
```

Because double precision allocates 52 bits to the significand, the default value of machine epsilon is $2^{-52}$.

```{code-cell}
log2(eps)
```

The spacing between adjacent floating-point values is proportional to the magnitude of the value itself. This is how relative precision is kept roughly constant throughout the range of values. You can get the adjusted spacing by calling `eps` with a value.

```{code-cell}
eps(1.618)
```

```{code-cell}
eps(161.8)
```

```{code-cell}
x = 161.8 + 0.1*eps(161.8);
x - 161.8
```

A common mistake is to think that $\epsilon_\text{mach}$ is the smallest floating-point number. It's only the smallest *relative to 1*. The correct perspective is that the scaling of values is limited by the exponent, not the mantissa. The actual range of positive values in double precision is

```{code-cell}
format short e
[realmin, realmax]
```

