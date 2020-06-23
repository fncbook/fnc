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

# Absolute and relative accuracy

Recall the grade-school approximation to the number $\pi$.

```{code-cell}
@show p = 22/7;
```

Not all the digits displayed for `p` are the same as those of $\pi$. As an approximation, its absolute and relative accuracy are

```{code-cell}
@show abs_accuracy = abs(p-pi);
```

```{code-cell}
@show rel_accuracy = abs(p-pi)/pi;
```

Note that $\pi$ is predefined in Julia. (You can also use the Greek letter for it by typing `\pi` and then pressing Tab.) Here we calculate the number of accurate digits in `p`.

```{code-cell}
@show accurate_digits = -log(10,rel_accuracy);
```
