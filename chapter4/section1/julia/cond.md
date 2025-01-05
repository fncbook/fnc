---
kernelspec:
  display_name: Julia 1
  language: julia
  name: julia-1.11
numbering:
  headings: false
---
```{code-cell}
include("../../../julia/FNC_init.jl")
```
[**Demo %s**](#demo-roots-cond)

Consider first the function

```{code-cell}
f(x) = (x - 1) * (x - 2);
```

:::{index} ! Julia; splatting
:::


At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at every point by a small amount of noise, we can imagine finding the root of the function drawn with a thick ribbon, giving a range of potential roots.
```{tip}
:class: dropdown
The syntax `interval...` is called **splatting** and means to insert all the individual elements of `interval` as a sequence.
```


```{code-cell}
interval = [0.8, 1.2]

plot(f, interval..., ribbon=0.03, aspect_ratio=1,
    xlabel=L"x", yaxis=(L"f(x)", [-0.2, 0.2]))

scatter!([1], [0], title="Well-conditioned root")
```

The possible values for a perturbed root all lie within the interval where the ribbon intersects the $x$-axis. The width of that zone is about the same as the vertical thickness of the ribbon.

By contrast, consider the function

```{code-cell}
f(x) = (x - 1) * (x - 1.01);
```

Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
plot(f, interval..., ribbon=0.03, aspect_ratio=1,
    xlabel=L"x", yaxis=(L"f(x)", [-0.2, 0.2]))

scatter!([1], [0], title="Poorly-conditioned root")
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function entirely upward by the amount drawn here, the root disappears!

