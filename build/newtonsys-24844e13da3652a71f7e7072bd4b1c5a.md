---
numbering:
  enumerator: 4.5.%s
---
(section-nonlineqn-newtonsys)=
# Newton for nonlinear systems

The rootfinding problem becomes much more difficult when multiple variables and equations are involved. 

```{index} ! rootfinding problem; multidimensional
```

````{prf:definition} Multidimensional rootfinding problem
Given a continuous vector-valued function $\mathbf{f}$ mapping from $\mathbb{R}^n$ into $\mathbb{R}^n$, find a vector $\mathbf{r}$ such that

```{math}
:label: rootvector
\begin{split}
  f_1(r_1,\dots,r_n) &= 0,\\
  f_2(r_1,\dots,r_n) &= 0,\\
  &\vdots\\
  f_n(r_1,\dots,r_n) &= 0.
\end{split}
```
````

Particular problems are often posed using scalar variables and equations. 

(example-newtonsys-predprey)=
::::{prf:example}
The steady state of interactions between the population $w(t)$ of a predator species and the population $h(t)$ of a prey species might be modeled as

$$
ah - b h w &= 0, \\ 
-cw + d w h &= 0
$$

for positive parameters $a,b,c,d$. To cast this in the form of {eq}`rootvector`, we could define $\mathbf{x}=[h,w]$, $f_1(x_1,x_2) = ax_1 - bx_1x_2$, and $f_2(x_1,x_2)= -c x_2 + d x_1 x_2$.
::::

While the equations of {numref}`Example {number} <example-newtonsys-predprey>` are easy to solve by hand, in practice even establishing the existence and uniqueness of solutions for any particular system is typically quite difficult.

## Linear model

To extend rootfinding methods to systems, we will keep to the basic philosophy of constructing easily managed models of the exact function. As usual, the starting point is a linear model. We base it on the multidimensional Taylor series,

```{math}
:label: multitaylor
\mathbf{f}(\mathbf{x}+\mathbf{h}) = \mathbf{f}(\mathbf{x}) + \mathbf{J}(\mathbf{x})\mathbf{h} + O(\| \mathbf{h} \|^2),
```

```{index} ! Jacobian matrix
```

where $\mathbf{J}$ is called the **Jacobian matrix** of $\mathbf{f}$ and is defined by

```{math}
:label: jacobian
\mathbf{J}(\mathbf{x}) =
  \begin{bmatrix}
    \rule[2mm]{0pt}{1em}\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\[2mm]
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\[1mm]
    \vdots & \vdots & & \vdots\\[1mm]
    \rule[-3mm]{0pt}{1em} \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
  \end{bmatrix} = \left[ \frac{\partial f_i}{\partial x_j} \right]_{\,i,j=1,\ldots,n}.
```

Because of the Jacobian's role in {eq}`multitaylor`, we may write $\mathbf{J}(\mathbf{x})$ as $\mathbf{f}{\,}'(\mathbf{x})$. Like any derivative, it is a function of the independent variable $\mathbf{x}$.

(example-nonlinsystem)=

````{prf:example}
Let
  
```{math}
\begin{split}
    f_1(x_1,x_2,x_3) &= -x_1\cos(x_2) - 1\\
    f_2(x_1,x_2,x_3) &= x_1x_2 + x_3\\
    f_3(x_1,x_2,x_3) &= e^{-x_3}\sin(x_1+x_2) + x_1^2 - x_2^2.
\end{split}
```

Then
  
```{math}
    \mathbf{J}(x) =
    \begin{bmatrix}
       -\cos(x_2) & x_1 \sin(x_2) & 0\\
      x_2 & x_1 & 1\\
       e^{-x_3}\cos(x_1+x_2)+2x_1 & e^{-x_3}\cos(x_1+x_2)-2x_2 &
       -e^{-x_3}\sin(x_1+x_2)
    \end{bmatrix}.
```

If we were to start writing out the terms in {eq}`multitaylor`, we would begin with
  
```{math}
\begin{split}
    f_1(x_1+h_1,x_2+h_2,x_3+h_3) &= -x_1\cos(x_2)-1 -\cos(x_2)h_1 +
    x_1\sin(x_2)h_2 + O\bigl(\| \mathbf{h} \|^2\bigr) \\
    f_2(x_1+h_1,x_2+h_2,x_3+h_3) &= x_1x_2 + x_3 + x_2h_1 +x_1h_2 +
    h_3 + O\bigl(\| \mathbf{h} \|^2\bigr),
  \end{split}
```

and so on.
````

The terms $\mathbf{f}(\mathbf{x})+\mathbf{J}(\mathbf{x})\mathbf{h}$ in {eq}`multitaylor` represent the linear part of $\mathbf{f}$ near $\mathbf{x}$. If $\mathbf{f}$ is actually linear, i.e., $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}-\mathbf{b}$, then the Jacobian matrix is the constant matrix $\mathbf{A}$ and the higher-order terms in {eq}`multitaylor` disappear.

## The multidimensional Newton iteration

```{index} Newton's method
```

With a method in hand for constructing a linear model for the vector system $\mathbf{f}(\mathbf{x})$, we can generalize Newton's method. Specifically, at a root estimate $\mathbf{x}_k$, we set $\mathbf{h} = \mathbf{x}-\mathbf{x}_k$ in {eq}`multitaylor` and get

```{math}
\mathbf{f}(\mathbf{x}) \approx \mathbf{q}(\mathbf{x})  = \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)(\mathbf{x}-\mathbf{x}_k).
```

We define the next iteration value $\mathbf{x}_{k+1}$ by requiring $\mathbf{q}(\mathbf{x}_{k+1})=\boldsymbol{0}$,

```{math}
\begin{split}
  \boldsymbol{0} &=  \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)(\mathbf{x}_{k+1}-\mathbf{x}_k),\\
\end{split}
```

which can be rearranged into

```{math}
:label: newtonsys
\mathbf{x}_{k+1} = \mathbf{x}_k - \bigl[\mathbf{J}(\mathbf{x}_k)\bigr]^{-1} \mathbf{f}(\mathbf{x}_k).
```

Note that $\mathbf{J}^{-1}\mathbf{f}$ now plays the role that $f/f'$ had in the scalar case; in fact, the two are the same in one dimension. In computational practice, however, we don't compute matrix inverses. 

```{index} ! Newton's method; multidimensional
```

(algorithm-nonlineqn-newtonsys)=
::::{prf:algorithm} Multidimensional Newton's method
Given $\mathbf{f}$ and a starting value $\mathbf{x}_1$, for each $k=1,2,3,\ldots$

1. Compute $\mathbf{y}_k = \mathbf{f}(\mathbf{x}_k)$ and $\mathbf{A}_k=\mathbf{f\,}'(\mathbf{x}_k)$.
2. Solve the linear system $\mathbf{A}_k\mathbf{s}_k = -\mathbf{y}_k$ for the **Newton step** $\mathbf{s}_k$.
3. Let $\mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{s}_k$.
::::

An extension of our series analysis of the scalar Newton's method shows that the vector version is also quadratically convergent in any vector norm, under suitable circumstances and when the iteration converges at all.

## Implementation

```{index} Julia; \\
```

An implementation of Newton's method for systems is given in {numref}`Function {number} <function-newtonsys>`. Other than computing the Newton step using backslash and taking vector magnitudes with `norm`, {numref}`Function {number} <function-newtonsys>` is virtually identical to the scalar version {numref}`Function {number} <function-newton>` presented earlier.

(function-newtonsys)=
``````{prf:algorithm} newtonsys
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-newtonsys-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-newtonsys-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-newtonsys-python
:::
````
`````
``````

(demo-newtonsys-converge)=
::::{prf:example}
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #demo-newtonsys-converge-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-newtonsys-converge-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-newtonsys-converge-python
:::
```` 
`````
::::

## Exercises

(problem-newtonsys-byhand)=
1. ✍ Suppose that
  
    ```{math}
    \mathbf{f}(\mathbf{x}) =
    \begin{bmatrix}
      x_1x_2+x_2^2-1 \\[1mm] x_1x_2^3 + x_1^2x_2^2 + 1
    \end{bmatrix}.
    ```

    Let $\mathbf{x}_1=[-2,1]^T$. Use Newton's method to find $\mathbf{x}_2$.

2. ✍ Suppose that $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x} - \mathbf{b}$ for a constant $n\times n$ matrix $\mathbf{A}$ and constant $n\times 1$ vector $\mathbf{b}$. Show that Newton's method converges to the exact root in one iteration.

    (problem-newtonsys-spherepotential)=
3. Two curves in the $(u,v)$ plane are defined implicitly by the equations $u\log u + v \log v = -0.3$ and $u^4 + v^2 = 1$.
  
    **(a)** ✍ Write the intersection of these curves in the form $\mathbf{f}(\mathbf{x}) = \boldsymbol{0}$ for two-dimensional $\mathbf{f}$ and $\mathbf{x}$.

    **(b)** ✍ Find the Jacobian matrix of $\mathbf{f}$.

    **(c)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find an intersection point starting from $u=1$, $v=0.1$.

    **(d)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find an intersection point starting from $u=0.1$, $v=1$.

    (problem-newtonsys-orbitintersect)=
4. Two elliptical orbits $(x_1(t),y_1(t))$ and $(x_2(t),y_2(t))$ are described by the equations
  
    ```{math}
    \begin{bmatrix}
      x_1(t) \\ y_1(t)
    \end{bmatrix}
    =
    \begin{bmatrix}
      -5+10\cos(t) \\ 6\sin(t)
    \end{bmatrix}, \qquad
    \begin{bmatrix}
      x_2(t)\\y_2(t)
    \end{bmatrix} =
    \begin{bmatrix}
      8\cos(t) \\ 3+12\sin(t)
    \end{bmatrix},
    ```

    where $t$ represents time.

    **(a)** ⌨ Make a plot of the two orbits with the following code:

    ``` julia
    x1(t) = -5+10*cos(t);   y1(t) = 6*sin(t);
    plot(x1,y1,0,2pi,aspect_ratio=1,legend=false)
    x2(t) = 8*cos(t);   y2(t) = 3 + 12*sin(t);
    plot!(x2,y2,0,2pi)
    ```

    **(b)** ✍ Write out a $2\times 2$ nonlinear system of equations that describes an intersection of these orbits. (Note: An intersection is not the same as a collision—they don't have to occupy the same point at the same time.)

    **(c)** ✍ Write out the Jacobian matrix of this nonlinear system.

    **(d)** ⌨ Use {numref}`Function {number} <function-newtonsys>` to find all of the unique intersections.
  
    (problem-newtonsys-ellipsemin)=
5. ⌨  Suppose one wants to find the points on the ellipsoid $x^2/25 + y^2/16 + z^2/9 = 1$ that are closest to and farthest from the point $(5,4,3)$. The method of Lagrange multipliers implies that any such point satisfies
  
    ```{math}
    \begin{split}
        x-5 &= \frac{\lambda x}{25}, \\[1mm]
        y-4 &= \frac{\lambda y}{16}, \\[1mm]
        z-3 &= \frac{\lambda z}{9}, \\[1mm]
        1 &=  \frac{1}{25}x^2 + \frac{1}{16}y^2 + \frac{1}{9}z^2
    \end{split}
    ```

    for an unknown value of $\lambda$.
  
    **(a)** Write out this system in the form $\mathbf{f}(\mathbf{u}) = \boldsymbol{0}$. (Note that the system has four variables to go with the four equations.)

    **(b)** Write out the Jacobian matrix of this system.

    **(c)** Use {numref}`Function {number} <function-newtonsys>` with different initial guesses to find the two roots of this system. Which is the closest point to $(5,4,3)$, and which is the farthest?
  
    (problem-newtonsys-circlefit)=
6. ⌨  Any three noncollinear points in the plane determine a unique circle. Suppose the points are given as $(x_i,y_i)$ for $i=1,2,3$. We can define the circle in terms of its center $(a,b)$ and radius $r$. Then 
    
    $$f_i(a,b,r) = (a-x_i)^2 + (b-y_i)^2 - r^2$$ 
    
    should be made zero for all $i=1,2,3$. This defines a nonlinear system $\mathbf{f}(\mathbf{v})=\boldsymbol{0}$ for $\mathbf{v}=[a,b,r]$. 

    Use {numref}`Function {number} <function-newtonsys>` on this system to find the circle passing through $(-5,0)$, $(1,-3)$, and $(4,2)$. Make a plot that shows you found the correct circle.



