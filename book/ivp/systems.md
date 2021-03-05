# Systems of differential equations

Before we improve on Euler's method, we will address a straightforward but critically important generalization of the problem. Very few applications involve an initial-value problem with just a single dependent variable. Most of the time there are multiple unknowns and a system of equations to define them.

````{prf:example}
Variations of the following model are commonly seen in ecology and epidemiology:

```{math}
  :label: predprey
  \begin{split}
    \frac{d y}{d t} &= y(1-\alpha y) - \frac{yz}{1+\beta y} \\
    \frac{d z}{d t} &= -z + \frac{yz}{1+\beta y},
  \end{split}
```

```{index} predator--prey model
```

where $\alpha$ and $\beta$ are positive constants. This model is a system of two differential equations for the unknown functions $y(t)$, which could represent a prey species or susceptible host, and $z(t)$, which could represent a predator species or infected population.  We refer to this as a **predator--prey model**. Both of the equations involve both of the unknowns, with no clear way to separate them.

We can pack the two dependent variables $y$ and $z$ into a vector-valued function of time, $\mathbf{u}(t)$, writing

```{math}
\begin{split}
  u_1'(t) &= f_1(t,\mathbf{u}) =  u_1(1-au_1) - \frac{u_1 u_2}{1+bu_1}\\
  u_2'(t) &= f_2(t,\mathbf{u}) = -u_2 + \frac{u_1 u_2}{1+bu_1}
\end{split}
```

and identifying $u_1=y$, $u_2=z$.
````

The generic form of a first-order system IVP is

```{math}
  :label: IVPsys
  \mathbf{u}'(t) = \mathbf{f}\bigl(t,\mathbf{u}(t)\bigr), \qquad a \le t \le b, \qquad
  \mathbf{u}(a)=\mathbf{u}_0,
```

```{prf:example} Julia demo
:class: demo
{doc}`demos/systems-exp`
```

where all of the boldface quantities are vectors in $\mathbb{R}^m$ or $\mathbb{C}^m$. This is an initial-value problem for a **first-order system of ODEs**, and the number $m$ of dependent variables (and equations) is called the **dimension** of the system. In particular, note that $\mathbf{f}$ represents $m$ scalar functions in $m+1$ scalar variables, including time.

```{prf:example} Julia demo
:class: demo
{doc}`demos/systems-predator`
```

The `DifferentialEquations` solvers all handle first-order systems. You must provide them with a function (either named or anonymous) that accepts the two arguments $t$ and $\mathbf{u}$ and returns $\mathbf{f}(t,\mathbf{u})$.

## Transformation of high-order systems

```{margin}
Every higher-order ODE system can be systematically transformed into a first-order system in a higher dimension.
```

Fortunately the ability to solve first-order ODE systems implies the ability to solve (most practical) systems of higher differential order, too. The reason is that there is a systematic way to turn a higher-order problem into a first-order one of higher dimension.

(example-nlosc3)=

````{prf:example}
 Consider the nonlinear initial-value problem
  
```{math}
  y''+(1+y')^3 y = 0, \qquad y(0)= y_0, \quad y'(0) = 0.
```

In order to write this problem as a first order system we define two scalar unknown functions, $u_1 = y$ and $u_2 = y'$. With these definitions, we have the two differential equations

```{math}
\begin{split}
  u_1' &= u_2 \\
  u_2' &= -(1+u_2)^3 u_1,
\end{split}
```

which is a first order system in two dimensions. The initial
condition of the system is

```{math}
  u_1(0) = y_0, \quad u_2(0) = 0.
```
````

(example-coupledpendula)=

````{prf:example}
Two identical pendula suspended from the same rod and swinging in parallel planes can be modeled as the second-order system

```{math}
\begin{split}
  \theta_1''(t) +\gamma \theta_1' + \frac{g}{L} \sin \theta_1 +
  k(\theta_1-\theta_2) &= 0\\
  \theta_2''(t) +\gamma \theta_2' + \frac{g}{L} \sin \theta_2 +
  k(\theta_2-\theta_1) &= 0,
\end{split}
```

where $\theta_1$ and $\theta_2$ are angles made by the two pendula, $L$ is the length of each pendulum, $\gamma$ is a frictional parameter, and $k$ is a parameter describing a torque produced by the rod when it is twisted. We can convert this problem into a first-order system using the substitutions

```{math}
  u_1 = \theta_1, \quad u_2 = \theta_2, \quad u_3 = \theta_1', \quad
  u_4 = \theta_2'.
```

With these definitions the system becomes

```{math}
\begin{split}
  u_1' &= u_3 \\
  u_2' &= u_4 \\
  u_3' &= -\gamma u_3 - \frac{g}{L}\sin u_1 + k(u_2-u_1) \\
  u_4' &= -\gamma u_4 - \frac{g}{L}\sin u_2 + k(u_1-u_2),
\end{split}
```

which is a first-order system in four dimensions. To complete the description of the problem, you would need to specify values for $\theta_1(0)$, $\theta_1'(0)$, $\theta_2(0)$, and $\theta_2'(0)$.
````

The trick illustrated in the preceding examples is always available. Specifically, one introduces a new variable (that is, a component of $\mathbf{u}$) for all but the highest derivative appearing for every variable of the original formulation. The surest way to get the transformation correct is to define one component of the new vector variable for each scalar initial condition given. Many equations for the first-order system then come from the trivial relationships among all the lower derivatives. The remaining equations for the system come from the original, high-order equations.  In the end, there must be as many scalar component equations as unknown first-order variables.

## Methods for IVP systems

```{index} Euler's method
```

The generalization of a scalar IVP solver to handle systems is straightforward. Consider Euler's method, which in system form becomes

```{math}
  :label: eulersys
  \begin{split}
    \mathbf{u}_{i+1} &= \mathbf{u}_i + h\mathbf{f}(t_i,\mathbf{u}_i), \qquad i=0,\ldots,n-1.
  \end{split}
```

The vector difference equation {eq}`eulersys` is just Euler's formula applied simultaneously to each component of the ODE system.  The method is still explicit for the solution at the new time level. Note here that as always in this book, an indexed boldface quantity such as $\mathbf{u}_i$ is a vector. If we want to refer to component $j$ of that vector, we would write $u_{i,j}$.

It might surprise you that the function {ref}`function-euler` that we wrote for scalar IVPs works for systems as well. In large part that's because as we observed with the mathematics, operations such as addition and multiplication translate easily from scalars to vectors. As long as the definition of the IVP correctly handles a vector for the dependent variable and specifies a vector initial condition, we're all set.

### new demo needed here

In the rest of this chapter we present methods as though they are for scalar equations, but their application to systems is straightforward. The generalization of error analysis can be more complicated, but our statements about order of accuracy and other properties are true for systems as well as scalars.  Our codes, on the other hand, are written to accept systems.

## Exercises

1. ✍ Rewrite the given higher order problems as first order systems.

    **(a)** $y'''-3y''+3 y' -y = t, \: y(0) = 1, \: y'(0) = 2, \: y''(0) = 3$

    **(b)** $y'' + 4 (x^2-1)y' + y = 0, \: y(0) = 2, \: y'(0) = -1$

    **(c)** For a given constant $a$,

    ```{math}
    \begin{split}
      x'' + \frac{a x}{(x^2+y^2)^{3/2}} &= 0\\
      y'' + \frac{a y}{(x^2+y^2)^{3/2}} &= 0,
      \end{split}
    ```

    with initial values $x(0) = 1$, $x'(0)=y(0) = 0$, $y'(0)=3$

    **(d)** $y^{(4)} -y = e^{-t}, \: y(0) = 0, \: y'(0) = 0, \: y''(0) = 1,\: y'''(0) = 0$

    **(e)** $y'''-y''+y'-y = t, \: y(0) = 1, \: y'(0) = 2, \: y''(0) = 3$

    ````{only} solutions
    ````

2. ✍ Write the given IVP as a system. Then do two steps of Euler's method by hand (perhaps with a calculator) with the indicated step size $h$. Using the given exact solution, compute the error after the second step.

    **(a)** $y''+ 4y = 4t, \: y(0) = 1,\: y'(0) = 1; \: \hat{y}(t) = t+\cos (2t),\: h=0.1$

    **(b)** $y''- 4y = 4t, \: y(0) = 2,\: y'(0) = -1; \: \hat{y}(t) = e^{2t} + e^{-2t}-t,\: h=0.1$

    **(c)** $2 x^2 y'' +3xy' - y = 0, \: y(2) = 1, \: y'(2) = -1/2,  \: \hat{y}(x) = 2/x, h = 1/8$

    **(d)** $2 x^2 y'' +3xy' - y = 0,\: y(1) = 4, \: y'(1) = -1, \: \hat{y}(x) = 2(x^{1/2} + x^{-1}), h=1/4$

    ````{only} solutions
    ````

3. ⌨ Solve the following IVPs using {ref}`function-euler` using $n=100$ steps. Plot the solution and its first derivative together on one plot, and plot the error in each component as functions of time on another.

    **(a)** $u''+ 4u = 4t, \: 0< t< 2\pi, \: u(0) = 1,\: u'(0) = 1; \: \hat{u}(t) = t+\cos (2t)$

    **(b)** $u''+ 9u = \sin(2t), \: 0< t< 2\pi, \: u(0) = 2,\: u'(0) = 1$;

    $\quad \hat{u}(t) = (1/5) \sin(3t) + 2 \cos (3t)+  (1/5) \sin (2t)$

    **(c)** $u''- 4u = 4t \: 0< t< \pi, \: u(0) = 2,\: u'(0) = -1; \: \hat{u}(t) = e^{2t} + e^{-2t}-t$

    **(d)** $u''+ 4u'+ 4u = t, \: 0< t< 4, \: u(0) = 1,\: u'(0) = 3/4; \: \hat{u}(t) = (3t+5/4)e^{-2t} + (t-1)/4$

    **(e)** $x^2 y'' +5xy' + 4y = 0,\: 1<x<e^2, \: y(1) = 0, \: y'(1) = 2, \: \hat{y}(x) = (2/x^2) \ln x$

    **(f)** $x^2 y'' +5xy' + 4y = 0,\: 1<x<e^2, \: y(1) = 1, \: y'(1) = -1, \: \hat{y}(x) = x^{-2}( 1 + \ln x)$

    **(g)** $2 x^2 y'' +3xy' - y = 0,\: 2<x<20, \: y(2) = 1, \: y'(2) = -1/2, \: \hat{y}(x) = 2/x$

    **(h)** $2 x^2 y'' +3xy' - y = 0,\: 1<x<16, \: y(1) = 4, \: y'(1) = -1, \: \hat{y}(x) = 2(x^{1/2} + x^{-1})$

    **(i)** $x^2 y'' -xy' + 2y = 0,\: 1<x<e^{\pi}, \: y(1) = 3, \: y'(1) = 4$;

    $\quad \hat{y}(x) = x \left[ 3 \cos \left( \ln x \right)+\sin \left( \ln x \right) \right]$

    **(j)** $x^2 y'' + 3xy' + 4y = 0,\: e^{\pi/12} < x < e^{\pi}, \: y(e^{\pi/12}) = 0,  \: y'(e^{\pi/12}) = -6$;

    $\quad \hat{y}(x) = x^{-1} \left[ 3 \cos \left( 3 \ln x \right)+\sin \left( 3 \ln x \right) \right]$

    ````{only} solutions
    ````

    (problem-SIR)=
4. ⌨ A disease that is endemic to a population can be modeled by tracking the fraction of the population that is susceptible to infection, $v(t)$, and the fraction that is infectious, $w(t)$. (The rest of the population is considered to be recovered and immune.) A typical model is the **SIR model** (see {cite}`brittonEssentialMathematical2003`)

    ```{math}
    \frac{dv}{dt} = 0.2(1-v) - 3vw, \qquad \frac{dw}{dt} = (3v-1)w.
    ```

    Starting with $v(0) = 0.95$ and $w(0) = 0.05$, use "ode45" to find the long-term steady values of $v(t)$ and $w(t)$. Plot both components of the solution as functions of time.

    ````{only} solutions
    f = @(t,u) [ 0.2*(1-u(1))-3*u(1)*u(2); 3*u(1)*u(2)-u(2) ];
    [t,u] = ode45(f,[0 50],[0.95;0.05]);
    plot(t,u)
    u(end,:)  % [ 0.3333 0.1333 ]
    ````

5. ⌨ Compute the solution to the systems for the given initial conditions using `solve`. Plot your results as a curve in the phase plane (that is, with $x$ and $y$ as the axes of the plot.)

6. Using initial conditions with $x(0)^2+y(0)^2$ both smaller and larger than 1 (inside and outside the unit circle), solve

    ```{math}
    \begin{split}
      x'(t) & = - 4y + x(1-x^2-y^2),\\
      y'(t) & = 4x + y(1-x^2-y^2),
    \end{split}
    ```

    starting at $0<t<10$.  What is the final state of the system?

7. Using initial conditions with $x(0)^2+y(0)^2$ both inside and outside circles of radius 1 and 2, solve

    ```{math}
    \begin{split}
      x'(t) & = - 4y + x(1-x^2-y^2)(4-x^2-y^2),\\
      y'(t) & = 4x + y(1-x^2-y^2)(4-x^2-y^2),
    \end{split}
    ```

    starting at $0<t<10$.  What is the final state of the system? Justify your answer with a plot showing the trajectories in the $(x,y)$-plane for a few different initial conditions.

    ````{only} solutions
    ````

    (problem-fitznag)=
8. ⌨ The **Fitzhugh--Nagumo equations** are a simple model of the repeated firing of a neuron. They are given by

    ```{math}
    \begin{split}
    \frac{d v_1}{dt} &= - v_1(v_1-1)(v_1-a) - v_2 + I, \\
    \frac{d v_2}{dt} &= \epsilon ( v_1 - \gamma v_2).
    \end{split}
    ```

    Assume $v_1(0) = 0.5$, $v_2(0) = 0.1$, $a = 0.1$, $\epsilon = 0.008$, $\gamma = 1$. For each value of $I$ below, find and plot the solution using `solve` for $0\le t \le 600$. The solutions are highly sensitive to $I$, and you need to change the requested absolute and relative error tolerances to $10^{-9}$. In each case the solution quickly approaches a periodic oscillation.

    **(a)** $I = 0.05527,\quad$
    **(b)** $I = 0.05683,\quad$
    **(c)** $I = 0.0568385,\quad$
    **(d)** $I = 0.05740$.

    This exploration was carried out by Baer and Erneux {cite}`baerSingularHopf1986`.

    ````{only} solutions
    a = 0.1;
    epsilon = 0.008;
    gamma = 1;

    I_ = [0.05527; 0.05683; 0.0568385; 0.05740];

    opt = odeset('reltol',1e-9,'abstol',1e-9);
    clf
    for I = I_'
        f = @(t,v) [ -v(1)*(v(1)-1)*(v(1)-a) - v(2) + I; ...
            epsilon*(v(1)-gamma*v(2)) ];
        [t,v] = ode45(f,[0,600],[0.5;0.1],opt);
        plot(t,v(:,1)), grid on, hold on
    end

    %%
    sol = ode45(f,[0,600],[0.5;0.1],opt);
    t1 = fzero(@(t) deval(sol,t,1), 300 );
    t2 = fzero(@(t) deval(sol,t,1)-0.1, 450 );
    period = t2-t1
    ````
