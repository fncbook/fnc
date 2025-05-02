---
numbering:
  enumerator: 12.4.%s
---
(section-advection-wave)=
# The wave equation

A close cousin of the advection equation is the **wave equation**.

```{index} ! wave equation
```

````{prf:definition} Wave equation
:::{math}
:label: wavepde
    u_{tt} - c^2 u_{xx} = 0.
:::
:::
````

The wave equation describes the propagation of waves, such as for sound or light. This is our first PDE having a second derivative in time. Consequently, we should expect to have two initial conditions: 

:::{math}
:label: waveIC
\begin{split}
u(x,0) &= f(x), \qquad 0 \le x \le 1,  \\
u_t(x,0) &= g(x), \qquad 0 \le x \le 1. 
\end{split}
:::

There are also two space derivatives, so if boundaries are present, two boundary conditions are required.

In the absence of boundaries, $u(x,t)=\phi(x-ct)$ is a solution of {eq}`wavepde`, just as it is for the advection equation. Now, though, so is $u(x,t)=\phi(x+c t)$ for any twice-differentiable $\phi$ (see @problem-wave-twodir).

```{note}
The wave equation in one dimension supports advection in both directions simultaneously.
```

<!-- We will use $x \in [0,1]$ as the domain.  -->


<!-- Because $u$ has two derivatives in $t$ and in $x$, we need two boundary conditions. We will use the Dirichlet conditions  

:::{math}
:label: waveBC
u(0,t) = u(1,t) = 0, \qquad t \ge 0,
::: -->

## First-order system

In order to apply semidiscretization with the standard IVP solvers that we have encountered, we must recast {eq}`wavepde` as a first-order system in time. Using our typical methodology, we would define $y=u_t$ and derive

:::{math}
:label: wavefirst1
\begin{split}
  u_t &= y, \\
  y_t &= c^2 u_{xx}.
\end{split}
:::

This is a valid approach. However, the wave equation has another, less obvious option for transforming to a first-order system:

:::{math}
:label: wavefirst2
\begin{split}
    u_t &= z_x, \\
    z_t &= c^2 u_{x}.
\end{split}
:::

```{index} Maxwell's equations
```

This second form is appealing in part because it's equivalent to [Maxwell's equations](wiki:Maxwell%27s_equations) for electromagnetism. In this form, we typically replace the velocity initial condition in {eq}`waveIC` with a condition on $z$, 

:::{math}
:label: waveIC2
\begin{split}
u(x,0) &= f(x),  \\
z(x,0) &= g(x).
\end{split}
:::

This alternative to {eq}`waveIC` is useful in many electromagnetic applications, where $u$ and $z$ represent the electric and magnetic fields, respectively.

Because waves travel in both directions, there is no preferred upwind direction. This makes a centered finite difference in space appropriate. Before application of the boundary conditions, semidiscretization of {eq}`wavefirst2` leads to

:::{math}
:label: waveMOL
  \begin{bmatrix}
    \mathbf{u}'(t) \\[2mm]  \mathbf{z}'(t)
  \end{bmatrix}
  =
  \begin{bmatrix}
    \boldsymbol{0} & \mathbf{D}_x \\[2mm] c^2 \mathbf{D}_x & \boldsymbol{0}
  \end{bmatrix}
  \begin{bmatrix}
    \mathbf{u}(t) \\[2mm] \mathbf{z}(t)
  \end{bmatrix}.
:::

The boundary conditions {eq}`waveBC` suggest that we should remove both of the end values of $\mathbf{u}$ from the discretization, but retain all of the $\mathbf{z}$ values. We use $\mathbf{w}(t)$ to denote the vector of all the unknowns in the semidiscretization. When computing $\mathbf{w}'(t)$, we extract the $\mathbf{u}$ and $\mathbf{z}$ components, and we use dedicated functions for padding with the zero end values or chopping off the zeros as necessary.

(demo-wave-boundaries)=
::::{prf:example} Wave equation with boundaries

We solve the wave equation {eq}`wavefirst2` with speed $c=2$, subject to {eq}`waveBC` and initial conditions {eq}`waveIC2`.

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-wave-boundaries-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-wave-boundaries-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-wave-boundaries-python
:::
````
`````
::::

## Variable speed

An interesting situation is when the wave speed $c$ changes discontinuously, as when light passes from one material into another. For this we must replace the term $c^2$ in {eq}`waveMOL` with the matrix $\operatorname{diag}\bigl(c^2(x_0),\ldots,c^2(x_m)\bigr)$.

(demo-wave-speed)=
::::{prf:example} Wave equation with variable speed

We now use a wave speed that is discontinuous at $x=0$; to the left, $c=1$, and to the right, $c=2$. 

`````{tab-set}
````{tab-item} Julia
:sync: julia
:::{embed} #demo-wave-speed-julia
:::
````

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-wave-speed-matlab
:::
````

````{tab-item} Python
:sync: python
:::{embed} #demo-wave-speed-python
:::
````
`````

::::

## Exercises

``````{exercise}
:label: problem-wave-maxwell
✍ Consider the Maxwell equations {eq}`wavefirst2` with smooth solution $u(x,t)$ and $z(x,t)$.

**(a)** Show that $u_{tt} = c^2 u_{xx}$.

**(b)** Show that $z_{tt} = c^2 z_{xx}$.
``````

``````{exercise}
:label: problem-wave-twodir
✍ Suppose that $\phi(s)$ is any twice-differentiable function.

**(a)** Show that $u(x,t) = \phi(x-c t)$ is a solution of $u_{tt}=c^2u_{xx}$. (As in the advection equation, this is a traveling wave of velocity $c$.)

**(b)** Show that $u(x,t) = \phi(x+c t)$ is another solution of $u_{tt}=c^2u_{xx}$. (This is a traveling wave of velocity $-c$.)
``````

```{index} D'Alembert's solution
```

``````{exercise}
:label: problem-wave-dalembert
✍ Show that the following is a solution to the wave equation $u_{tt}=c^2u_{xx}$ with initial and boundary conditions {eq}`waveBC` and {eq}`waveIC`:

$$
u(x,t) = \frac{1}{2} \left[ f(x-ct)+f(x+ct)\right] + \frac{1}{2c} \int_{x-ct}^{x+ct} g(\xi) \, d\xi
$$

This is known as [D'Alembert's formula](wiki:D%27Alembert%27s_formula).
``````

``````{exercise}
:label: problem-wave-neumann
⌨ Suppose the wave equation has homogeneous Neumann conditions on $u$ at each boundary instead of Dirichlet conditions. Using the Maxwell formulation {eq}`wavefirst2`, we have $z_t=c^2u_x$, so $z$ is constant in time at each boundary. Therefore, the endpoint values of $\mathbf{z}$ can be taken from the initial condition and removed from the ODE, while the entire $\mathbf{u}$ vector is now part of the ODE. 

Modify @demo-wave-boundaries to solve the PDE there with Neumann instead of Dirichlet conditions, and make an animation or space-time portrait of the solution. In what major way is it different from the Dirichlet case?
``````

```{index} sine–Gordon equation
```

``````{exercise}
:label: problem-wave-sinegordon
The nonlinear [sine–Gordon equation](wiki:Sine-Gordon_equation) $u_{tt}-u_{xx}= - \sin u$ is a model of multiple pendulums hanging from a single string. It has some interesting solutions.

**(a)** ✍ Write the equation as a first-order system in the variables $u$ and $v=u_t$.

**(b)** ⌨ Assume periodic end conditions on $[-10,10]$ and discretize at $m=200$ points. Let 

$$
u(x,0) = 4\operatorname{arctan}\left( \frac{1}{\cosh(x/\sqrt{2})} \right) 
$$ 

and $u_t(x,0) = 0$. Solve the system using a nonstiff solver for $t \in [0,50],$ and make an animation of the solution $u(x,t)$. (This type of solution is called a [breather](wiki:Breather).)
``````

```{index} beam equation
```

``````{exercise}
:label: problem-wave-beam

The deflections of a stiff beam, such as a ruler, are governed by the PDE $u_{tt}=-u_{xxxx}$.

**(a)** ✍ Show that the beam PDE is equivalent to the first-order system

\begin{align*}
u_t &= v_{xx}, \\
v_t &= -u_{xx}.
\end{align*}

**(b)** ⌨ Assuming periodic end conditions on $[-1,1]$, use $m=100$, let $u(x,0) =\exp(-24x^2)$, $v(x,0) = 0$, and simulate the solution of the beam equation for $0\le t \le 1$ using {numref}`Function {number} <function-am2>` with $n=100$ time steps. Make a plot or animation of the solution.
``````
