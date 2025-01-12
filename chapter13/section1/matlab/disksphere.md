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
[**Demo %s**](#demo-tensorprod-disksphere)

For a function given in polar form, such as $f(r,\theta)=1-r^4$, construction of a function over the unit disk is straightforward using a grid in $(r,\theta)$ space.

```{code-cell}
r = linspace(0, 1, 41);
theta = linspace(0, 2*pi, 121);
[mtx, R, Theta] = tensorgrid(r, theta);
F = mtx(@(r, theta) 1 - r.^4);
clf,  colormap(parula)
contourf(R', Theta', F', 20)
shading flat
xlabel("r"),  ylabel("\theta"), 
title("A polar function")   
```

Of course, we are used to seeing such plots over the $(x,y)$ plane, not the $(r,\theta)$ plane. For this we create matrices for the coordinate functions $x$ and $y$.

```{code-cell}
X = R .* cos(Theta);  Y = R .* sin(Theta);
contourf(X', Y', F', 20)
axis equal,  shading interp  
xlabel('x'),  ylabel('y')
title('Function over the unit disk')  
```

In such functions the values along the line $r=0$ must be identical, and the values on the line $\theta=0$ should be identical to those on $\theta=2\pi$. Otherwise the interpretation of the domain as the unit disk is nonsensical. If the function is defined in terms of $x$ and $y$, then those can be defined in terms of $r$ and $\theta$ using {eq}`unitdiskparam`.

