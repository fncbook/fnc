---
numbering:
  enumerator: 10.2.%s
---
(section-bvp-shooting)=
# Shooting

One way to attack the TPBVP {eq}`tpbvp` is to adapt our IVP solving techniques from [Chapter 6](../ivp/overview.md) to it. Those techniques work only when we know the entire initial state, but we can allow that state to vary in order to achieve the stated conditions. 

This is the idea behind the **shooting method**. Imagine adjusting your aiming point and power to sink a basketball shot from the free-throw line. The way in which you miss—too long, flat trajectory, etc.—informs how you will adjust for your next attempt.

(demo-shooting-naive)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Naive shooting
:open:
```{include} julia/naive.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Naive shooting
:open:
```{include} matlab/naive.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Naive shooting
:open:
```{include} python/naive.ipynb
```
````
`````
``````
```````
    

We can do much better than trial-and-error for the unknown part of the initial state. As usual, we can rewrite the ODE $u''(x) = \phi(x,u,u')$ in first-order form as

:::{math}
:label: shoot-first
\begin{split}
y_1' &= y_2,\\ 
y_2' &= \phi(x,y_1,y_2).
\end{split}
:::

We turn this into an IVP by specifying $y(a)=s_1$, $y'(a)=s_2$, for a vector $\mathbf{s}$ to be determined by the boundary conditions. Define the residual function $\mathbf{v}(\mathbf{s})$ by

:::{math}
:label: shoot-resid
\begin{split}
v_1(s_1,s_2) &= g_1(y_1(a),y_2(a)) = g_1(s_1,s_2),\\ 
v_2(s_1,s_2) &= g_2(y_1(b),y_2(b)).
\end{split}
:::

The dependence of $v_2$ on $\mathbf{s}$ is indirect, through the solution of the IVP for $\mathbf{y}(x)$. We now have a standard rootfinding problem that can be solved via the methods of [Chapter 4](../nonlineqn/overview.md). 

## Implementation

Our implementation of shooting is given in {numref}`Function {number} <function-shoot>`. Note the structure: we use a rootfinding method that in turn relies on an IVP solver. This sort of arrangement is what makes us concerned with minimizing the number of objective function calls when rootfinding.

(function-shoot)=
``````{prf:algorithm} shoot
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{embed} #function-shoot-julia
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #function-shoot-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #function-shoot-python
:::
````
`````
``````

(demo-shooting-mems)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Shooting solution of a BVP
:open:
```{include} julia/mems.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Shooting solution of a BVP
:open:
```{include} matlab/mems.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Shooting solution of a BVP
:open:
```{include} python/mems.ipynb
```
````
`````
``````
```````
    

## Instability

The accuracy of the shooting method should be comparable to those of the component pieces, the rootfinder, and the IVP solver. However, the shooting method is unstable for some problems. An example illustrates the trouble.

(demo-shooting-unstable)=
```````{prf:example}
``````{tab-set}
`````{tab-item} Julia
:sync: julia
````{dropdown} Instability of shooting
:open:
```{include} julia/unstable.ipynb
```
````
`````

`````{tab-item} MATLAB
:sync: matlab
````{dropdown} Instability of shooting
:open:
```{include} matlab/unstable.ipynb
```
````
`````

`````{tab-item} Python
:sync: python
````{dropdown} Instability of shooting
:open:
```{include} python/unstable.ipynb
```
````
`````
``````
```````
    


The essence of the instability is that errors can grow exponentially away from the boundary at $x=a$, where the state is arbitrarily being set (see {numref}`Theorem {number} <theorem-depIC>`). Using shooting, acceptable accuracy near $x=b$ therefore means requiring extraordinarily high accuracy near $x=a$.

The instability of shooting can be circumvented by breaking the interval into smaller pieces and thus limiting the potential for error growth. However, we do not go into these details. Instead, the methods in the rest of this chapter treat both ends of the domain symmetrically and solve over the whole domain at once.

## Exercises

1. ⌨ For each BVP in [Exercise 10.1.2](#problem-tpbvp-verify) , use {numref}`Function {number} <function-shoot>` to compute the solution. Plot the solution and, separately, its error as functions of $x$. 

    ```{index} pendulum
    ```
(problem-shooting-pendulum)=
2. ⌨ (Continuation of [Exercise 10.1.4](#problem-tpbvp-allencahn).) Consider the pendulum from {numref}`Example {number} <example-tpbvp-pendulum>` with $g=L=1$. Suppose we want to release the pendulum from rest such that $\theta(5)=\pi/2$. Find one solution that passes through $\theta=0$, and another solution that does not. Plot $\theta(t)$ for both cases together.

    ```{index} Allen–Cahn equation
    ```
(problem-shooting-allencahn)=
3. ⌨  (Continuation of [Exercise 10.1.5](#problem-tpbvp-allencahn).) The stationary Allen–Cahn equation is 
 
    $$
      \epsilon u'' = u^3-u, \qquad 0 \le x \le 1, \qquad u(0)=-1, \quad u(1)=1.
    $$

    As $\epsilon\rightarrow 0$, the solution tends toward a step function transition between $-1$ and $1$. By symmetry, $u'(x)=-u'(1-x)$.
  
    **(a)** Use {numref}`Function {number} <function-shoot>` to solve the equation for $\epsilon=0.2$. Plot the solution and compute the numerical value of $u'(0)-u'(1)$.
    
    **(b)** Repeat for $\epsilon=0.02$.
    
    **(c)** Repeat for $\epsilon=0.002$. You will receive multiple warning messages. Does the result look like a valid solution?

4. ✍ Consider the linear TPBVP 
    
    $$
    \begin{split}
    u'' &= p(x)u' + q(x)u + r(x),\\ 
    u'(a) &= 0, \quad u(b)=\beta.
    \end{split}
    $$

    The shooting IVP uses the same ODE with initial data $u(a)=s_1$, $u'(a)=s_2$ to solve for a trial solution $u(x)$. Define

    $$
    z(x) = \frac{\partial u}{\partial s_1}.
    $$

    By differentiating the IVP with respect to $s_1$, show that $z$ satisfies the IVP

    $$
    z'' = p(x)z' + q(x)z, \quad z(0)=1, \; z'(0)=0.
    $$

    It follows that $z(x)$ is independent of $s_1$, and therefore $u(x)$ is a linear function of $s_1$ at each fixed $x$. Use the same type of argument to show that $u(x)$ is also a linear function of $s_2$, and explain why the residual function $\mathbf{v}$ in {eq}`shoot-resid` is a linear function of $\mathbf{s}$.
