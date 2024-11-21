---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---
# Advection equations

```{index} Obi-Wan Kenobi, Star Wars: A New Hope
```

```{epigraph}
Now, let's see if we can't figure out what you are, my little friend. And where you come from.

— Obi-Wan Kenobi, *Star Wars: A New Hope* 
```

Now that we have seen PDEs with both time and space variables, we have a new wrinkle to add. Some of these equations behave like those of the previous chapter, creating diffusive effects. Others, though, are about propagation or advection—generally, the behavior of waves.

Wave behavior is very different from diffusion. In idealized cases, waves travel with a finite speed and conserve energy, whereas diffusion smooths out features quickly and is associated with the dissipation of energy. There are many numerical methods that are specialized for purely advective problems, but we will consider only the method of lines approach from [Chapter 11](../diffusion/overview). Along the way we will see that wave-like behavior leads to some different conclusions about how to make these methods effective.

