# Newton for nonlinear systems

The rootfinding problem becomes much more difficult when multiple variables and equations are involved. We now let $\mathbf{f}$ be a vector function of a vector argument: $\mathbf{f}$ maps from $\mathbb{R}^n$ to $\mathbb{R}^n$. By $\mathbf{f}(\mathbf{x})=\boldsymbol{0}$ we mean the simultaneous system of $n$ scalar equations,

```{math}
\begin{split}
  f_1(x_1,\dots,x_n) &= 0\\
  f_2(x_1,\dots,x_n) &= 0\\
  &\vdots\\
  f_n(x_1,\dots,x_n) &= 0.
\end{split}
```

When discussing a specific problem, it is often necessary to write out the equations componentwise, but in discussing methods for the general problem it's more convenient to use the vector form. Proving the existence and uniqueness of a solution for any particular $\mathbf{f}$ is typically quite difficult.

## Linear model

To extend rootfinding methods to systems, we will keep to the basic philosophy of constructing easily managed models of the exact function. As usual, the starting point is a linear model. We base it on the multidimensional Taylor series,

```{math}
:label: multitaylor
\mathbf{f}(\mathbf{x}+\mathbf{h}) = \mathbf{f}(\mathbf{x}) + \mathbf{J}(\mathbf{x})\mathbf{h} + O(\| \mathbf{h} \|^2),
```

```{index} Jacobian matrix
```

where $\mathbf{J}$ is called the {term}`Jacobian matrix` of $\mathbf{f}$ and is defined by

```{math}
:label: jacobian
\mathbf{J}(\mathbf{x}) =
  \begin{bmatrix}
    \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n}\\[2mm]
    \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n}\\[1mm]
    \vdots & \vdots & & \vdots\\[1mm]
    \frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \cdots & \frac{\partial f_n}{\partial x_n}
  \end{bmatrix} = \left[ \frac{\partial f_i}{\partial x_j} \right]_{i,j=1,\ldots,n.}
```

Because of the Jacobian's role in {eq}`multitaylor`, we may write $\mathbf{J}(\mathbf{x})$ as $\mathbf{f}\,'(\mathbf{x})$. Like any derivative, it is a function of the independent variable $\mathbf{x}$.

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
    J(x) =
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

The terms $\mathbf{f}(\mathbf{x})+\mathbf{J}(\mathbf{x})\mathbf{h}$ in {eq}`multitaylor` represent the "linear part" of $\mathbf{f}$ near $\mathbf{x}$. If $\mathbf{f}$ is actually linear, i.e., $\mathbf{f}(\mathbf{x})=\mathbf{A}\mathbf{x}-\mathbf{b}$, then the Jacobian matrix is the constant matrix $\mathbf{A}$ and the higher order terms in {eq}`multitaylor` disappear.

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

Note that $\mathbf{J}^{-1}\mathbf{f}$ now plays the role that $f/f'$ had in the scalar case; in fact the two are the same in one dimension. In computational practice, however, we don't compute matrix inverses. Instead, define the *Newton step* $\mathbf{s}_k$ by the linear $n\times n$ system

```{math}
  :label: newtonstep
  \mathbf{J}(\mathbf{x}_k)\, \mathbf{s}_k = -\mathbf{f}(\mathbf{x}_k),
```

```{margin}
Computing the Newton step is equivalent to solving a linear system using the Jacobian matrix and the function value.
```

so that $\mathbf{x}_{k+1}=\mathbf{x}_k+\mathbf{s}_k$. Computing the Newton step is equivalent to solving a linear system using the Jacobian matrix and the function value.

````{prf:example} Julia demo
:class: demo
{doc}`demos/system-iter`
````

An extension of our series analysis of the scalar Newton's method shows that the vector version is also quadratically convergent in any vector norm, under suitable circumstances.

## Implementation

An implementation of Newton's method for systems is given in {ref}`function-newtonsys`.

(function-newtonsys)=

```{proof:function} newtonsys
**Newton's method for a system of equations.**

```{code-block} julia
:lineno-start: 1
"""
newtonsys(f,jac,x1)

Use Newton's method to find a root of a system of equations,
starting from `x1`. The functions `f` and `jac should return the
residual vector and the Jacobian matrix, respectively. Returns
history of root estimates as a vector of vectors.
"""
function newtonsys(f,jac,x1)
    # Operating parameters.
    funtol = 1000*eps();  xtol = 1000*eps();  maxiter = 40;

    x = [float(x1)]
    y,J = f(x1),jac(x1)
    dx = Inf   # for initial pass below
    k = 1

    while (norm(dx) > xtol) && (norm(y) > funtol) && (k < maxiter)
        dx = -(J\y)             # Newton step
        push!(x,x[k] + dx)    # append to history
        k += 1
        y,J = f(x[k]),jac(x[k])
    end

    if k==maxiter
        @warn "Maximum number of iterations reached."
    end

    return x
end
```

````{tip}
```{toggle}
The output of {ref}`function-newtonsys` is a vector of vectors representing the entire history of root estimates. Since these are usually floating point, even when the initial estimate is integer, the initial estimate is converted with `float` in line 12.
```
````

````{prf:example} Julia demo
:class: demo
{doc}`demos/system-usage`
````

It is a remarkable effect of the vector-friendliness of Julia that this program is hardly different from the scalar version {ref}`function-newton` presented earlier:

```{index} backslash
```

- The root estimates are stored as columns in an array.
- The Newton step is calculated using a backslash.
- The function "norm" is used for the magnitude of a vector, instead of "abs" for the magnitude of a scalar.

Indeed, {ref}`function-newtonsys` is a proper generalization—it can be used on scalar problems as well as on systems.

## Exercises

(problem-newtonsysbyhand)=

1. ✍  Suppose that
  
    ```{math}
    \mathbf{f}(\mathbf{x}) =
    \begin{bmatrix}
      x_1x_2+x_2^2-1 \\[1mm] x_1x_2^3 + x_1^2x_2^2 + 1
    \end{bmatrix}.
    ```

    Let $\mathbf{x}_1=[-2,1]^T$. Use Newton's method to find $\mathbf{x}_2$.

    ````{only} solutions
    %% (a) $f(x_1)=[-2;3]$ and $J(x_1)=[ 1,0; -3,2]$. The Newton step is $-[-2;-3/2]$ and the result is $[0;5/2]$.
    ````

2. ✍ Suppose that $\mathbf{f}(\mathbf{x}) = \mathbf{A}\mathbf{x} - \mathbf{b}$ for a constant $n\times n$ matrix $\mathbf{A}$ and constant $n\times 1$ vector $\mathbf{b}$. Show that Newton's method converges to the exact root in one iteration.

    (problem-spherepotential)=
3. Two curves in the $(u,v)$ plane are defined implicitly by the equations $u\log u + v \log v = -0.3$ and $u^4 + v^2 = 1$.
  
    **(a)** ✍ Write the intersection of these curves in the form $\mathbf{f}(\mathbf{x}) = \boldsymbol{0}$ for two-dimensional $\mathbf{f}$ and $\mathbf{x}$.

    **(b)** ✍ Find the Jacobian matrix of $\mathbf{f}$.

    **(c)** ⌨ Use {ref}`function-newtonsys` to find an intersection point near $u=1$, $v=0.1$.

    **(d)** ⌨ Use {ref}`function-newtonsys` to find an intersection point near $u=0.1$, $v=1$.

    ````{only} solutions

    ``` matlab
    f = @(x) deal( [x(1)*log(x(1))+x(2)*log(x(2))+0.3;x(1)^4+x(2)^2-1],...
        [log(x(1))+1,log(x(2))+1;4*x(1)^3,2*x(2)^1] );
    newtonsys(f,[1;.1])
    newtonsys(f,[.1;1])
    ```
    ````

    (problem-orbitintersect)=
4. Two elliptical orbits $(x_1(s),y_1(s))$ and $(x_2(t),y_2(t))$ are described by the equations
  
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
      8\cos(t) \\ 1+12\sin(t)
    \end{bmatrix},
    ```

    where $t$ represents time.

    **(a)** ⌨ Make a plot of the two orbits in the following code:

    ``` julia
    x1(t) = -5+10*cos(t);   y1(t) = 6*sin(t);
    plot(x1,y1,0,2pi,aspect_ratio=1,legend=false)
    x2(t) = 8*cos(t);   y2(t) = 1+12*sin(t);
    plot!(x2,y2,0,2pi)
    ```

    **(b)** ✍ Write out a $2\times 2$ nonlinear system of equations that describes an intersection of these orbits. (Note: An intersection is not the same as a collision—they don't have to occupy the same point at the same time.)

    **(c)** ✍ Write out the Jacobian matrix of this nonlinear system.

    **(d)** ⌨ Use {ref}`function-newtonsys` to find all of the unique intersections.
  
    ````{only} solutions
    ``` matlab
    clf
    x1 = @(t) -5+10*cos(t);   y1 = @(t) 6*sin(t);
    fplot(x1,y1,[0 2*pi])
    hold on
    x2 = @(t) 8*cos(t);   y2 = @(t) 1+12*sin(t);
    fplot(x2,y2,[0 2*pi])


    f = @(x) deal( [-5+10*cos(x(1)) - 8*cos(x(2)); 6*sin(x(1))-1-12*sin(x(2))],...
        [-10*sin(x(1)),8*sin(x(2)); 6*cos(x(1)),-12*cos(x(2))] );
    x = newtonsys(f,[-8;9]);
    t = x(:,end);
    pos1 = [x1(t(1));y1(t(1))]
    x = newtonsys(f,[8;-9]);
    t = x(:,end);
    pos2 = [x1(t(1));y1(t(1))]
    ```
    ````

    (problem-ellipsemin)=
5. ⌨  Suppose one wants to find the points on the ellipsoid $x^2/25 + y^2/16 + z^2/9 = 1$ that are closest to and farthest from the point $(5,4,3)$. The method of Lagrange multipliers implies that any such point satisfies
  
    ```{math}
    \begin{split}
        x-5 &= \frac{\lambda x}{25} \\
        y-4 &= \frac{\lambda y}{16} \\
        z-3 &= \frac{\lambda z}{9} \\
        1 &=  \frac{1}{25}x^2 + \frac{1}{16}y^2 + \frac{1}{9}z^2
    \end{split}
    ```

    for an unknown value of $\lambda$.
  
    **(a)** Write out this system in the form $\mathbf{f}(\mathbf{u}) = \boldsymbol{0}$. (Note that the system has four variables to go with the four equations.)

    **(b)** Write out the Jacobian matrix of this system.

    **(c)** Use {ref}`function-newtonsys` with different initial guesses to find the two roots of this system. Which is the closest point to $(5,4,3)$ and which is the farthest?
  
    ````{only} solutions
    ``` matlab
    %%
    %

    %% part(a)
    % Let $u_1=x$, $u_2=y$, $u_3=z$, $u_4=\lambda$. Then $f_1(\mathbf{u})=u_1-5-u_4 u_1 / 25$, etc.

    %% part(b)
    % $$J(\mathbf{u}) =[   1 - u_4/25, 0 ,0 , -u_1/25; 0, 1 - u_4/16 ,0 , -u_2/16;$$

    %%
    % $$0, 0, 1-u_4/9 , -u_3/9;  2u_1/25 , 2u_2/16 , 2u_3/9 , 0 , 0 ]$$

    %% part (c)
    f = @(u) [ u(1)-5-u(1)*u(4)/25; ...
              u(2)-4-u(2)*u(4)/16; ...
              u(3)-3-u(3)*u(4)/9; ...
              u(1)^2/25 + u(2)^2/16 + u(3)^2/9 - 1 ];

    J = @(u) [ 1-u(4)/25, 0 ,0 , -u(1)/25; ...
              0, 1-u(4)/16 ,0 , -u(2)/16; ...
              0, 0, 1-u(4)/9 , -u(3)/9;...
              2*u(1)/25 , 2*u(2)/16 , 2*u(3)/9 , 0 ];

    %%
    u0 = [1;2;3;4];
    u = newtonsys(f,J,u0);
    point1 = u(:,end)

    %%
    u = newtonsys(f,J,-point1);
    point2 = u(:,end)
    ```
    ````

6. ⌨  In this problem you are to fit a function of the form
  
    ```{math}
    P(t) = a_1 + a_2 e^{a_3 t}
    ```

    to a subset of US census data for the twentieth century:

    | Year | 1910 | 1930 | 1950 | 1970 | 1990 |
    |------|------|------|------|------|------|
    | Population | 92.0 | 122.8 | 150.7 | 205.0$ | 248.7 |
  
    **(a)** Determine the unknown parameters $a_1$, $a_2$, $a_3$ in $P$ by requiring that $P$ exactly reproduce the data in the years 1910, 1950, and 1990. This creates three nonlinear equations for $a_1$, $a_2$, and $a_3$ that may be solved using {ref}`function-newtonsys`.

    **(b)** To obtain convergence, rescale the data using the time variable $t = (\text{year}-1900)/100$ and divide the population numbers above by  $100$. Using your model $P(t)$, predict the result of the 2000 census, and compare it to the true figure of 284.1 million.
  
