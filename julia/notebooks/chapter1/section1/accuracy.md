---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
Recall the grade-school approximation to the number $\pi$.

```{code-cell}
@show p = 22/7;
```
Not all the digits displayed for `p` are the same as those of $\pi$. 
```{tip}
The value of `pi` is predefined and equivalent to `π`, which is entered by typing `\pi` followed immediately by the <kbd>Tab</kbd> key.
```

```{code-cell}
@show float(π);
```

```{index} ! Julia; string interpolation
```


The absolute and relative accuracies of the approximation are as follows.
```{tip}
A dollar sign `$` in a string substitutes (or *interpolates*) the named variable or expression into the string.
```

```{code-cell}
acc = abs(p-π)
println("absolute accuracy = $acc")
println("relative accuracy = $(acc/π)")
```

Here we calculate the number of accurate digits in `p`.
```{tip}
The `log` function is for the natural log. For other common bases, use `log10` or `log2`.
```

```{code-cell}
println("Number of accurate digits = $(-log10(acc/π))")
```
This last value could be rounded down by using `floor`.

