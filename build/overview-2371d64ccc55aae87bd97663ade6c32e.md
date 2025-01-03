(chapter-ivp)=
# Initial-value problems for ODEs

```{index} Han Solo, Star Wars: A New Hope
```

```{epigraph}
Without precise calculations we could fly right through a star or bounce too close to a supernova and that'd end your trip real quick, wouldn't it?

— Han Solo, *Star Wars: A New Hope*
```

Quantities that change continuously in time or space are often modeled by differential equations. When everything depends on just one independent variable, we call the model an **ordinary differential equation** (ODE).  Differential equations need supplemental conditions to define both the modeling situation and the theoretical solutions uniquely. The *initial-value problem* (IVP), in which all of the conditions are given at a single value of the independent variable, is the simplest situation. Often the independent variable in this case represents time.

Methods for IVPs usually start from the known initial value and iterate or "march" forward from there. There is a large number of them, owing in part to differences in accuracy, stability, and convenience. The most broadly important methods fall into one of two camps: *Runge–Kutta* and *linear multistep* formulas. Each type introduces its own complications, and we will consider them separately.

