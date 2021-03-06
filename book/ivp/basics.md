# Basics of IVPs

```{index} initial-value problem
```

The general form of a scalar, first-order {term}`initial-value problem` (IVP) is
  
```{math}
:label: IVP
\begin{split}
   u'(t) &= f(t,u(t)), \qquad a \le t \le b  \\
  u(a) &=u_0.
\end{split}
```

We call $t$ the **independent variable** and $u$ the **dependent variable**. When $t$ is in fact meant to be time, sometimes we write $\dot{u}$ (read "u-dot") instead of $u'$. A solution of an initial-value problem is a function $u(t)$ that makes both $u'(t)=f\bigl(t,u(t)\bigr)$ and $u(a)=u_0$ true equations.

````{prf:example}
  Suppose $u(t)$ is the size of a population at time $t$. We idealize by allowing $u$ to take any real (not just integer) value. If we assume a constant per capita birth rate (births per unit population per unit time), then
  
```{math}
\frac{d u}{d t} = ku, \qquad u(0)=u_0,
```

for some $k>0$. The solution is $u(t)=e^{kt}u_0$, which is exponential growth.

A more realistic model would cap the growth due to finite resources. Suppose the death rate is proportional to the size of the population, indicating competition. Then

```{math}
  :label: logistic
  \frac{d u}{d t} = ku - ru^2, \qquad u(0)=u_0.
```

```{index} logistic equation
```

This is the **logistic equation**. Although crude, it is still useful in population models.  The solution relevant for population models has the form
  
```{math}
  u(t) = \frac{k/r}{ 1 + \left( \frac{k}{r u_0} - 1 \right) e^{-k t} }.
```

For $k,r,u_0>0$, the solution smoothly varies from the initial population $u_0$ to a finite population, equal to $k/r$, that has been limited by competition.
````

If $u'=f(t,u)=g(t)+u h(t)$, the differential equation is **linear**. Linear problems can be solved in terms of standard integrals. Defining the **integrating factor** $\rho(t) = \exp\bigl[\int -h(t)\, dt \bigr]$, the solution is derived from

```{math}
  \rho(t) u(t) = u_0 + \int_a^t \rho(s) g(s) \, ds.
```

In many cases, however, the integrals cannot be done in closed form. When the differential equation is nonlinear, there is usually no analytic formula available for its solution.

An ODE may have higher derivatives of the unknown solution present. For example, a **second-order ordinary differential equation** is often given in the form $u''(t)=f\bigl(t,u,u'\bigr)$. A second-order IVP requires two conditions at the initial time in order to specify a solution completely. As we will see in {doc}`systems`, we are always able to reformulate higher-order IVPs in a first-order form, so we will deal with first-order problems exclusively. (There are, however, some numerical methods for second-order problems that are important to certain application areas.)

## Numerical solutions

```{prf:example} Julia demo
:class: demo
:label: demos-basics-first
{doc}`demos/basics-first`

:label: demos-basics-usage
{doc}`demos/basics-usage`
```

The `DifferentialEquations` package has numerous methods for solving ordinary differential equations, including initial-value problems. Underlying each numerical solution is a set of point values of the solution. The package includes interpolation procedures that allow you to evaluate the solution at any time, however.

## Existence and uniqueness

```{prf:example} Julia demo
:class: demo
:label: demos-basics-sing
{doc}`demos/basics-sing`
```

As demonstrated in {doc}`demos/basics-sing`, there are simple IVPs that do not have solutions at all possible times. Furthermore, we can easily find an IVP that has more than one solution.

````{prf:example}
  The functions $u(t)=u^2$ and $u(t)\equiv 0$ both satisfy the differential equation $u'=2\sqrt{u}$ and the initial condition $u(0)=0$. Thus the corresponding IVP has more than one solution.
````

The following standard theorem gives us a condition that is easy to check and guarantees that a unique solution exists. But it is not the most general possible such condition, so there are problems with a unique solution that it cannot detect. We state the theorem without proof.

(theorem-existunique)=

````{prf:theorem} (Existence and uniqueness)
 If the derivative $\frac{\partial f}{\partial u}$ exists and $\left|\frac{\partial f}{\partial u}\right|$ is bounded by a constant $L$ for all $a\le t \le b$ and all $u$, then the initial-value problem {eq}`IVP` has a unique solution for $t\in [a,b]$.
````

## Conditioning of first-order IVPs

```{index} condition number; of initial-value problems
```

In a numerical context we have to be concerned about the conditioning of the IVP. There are two key items in {eq}`IVP` that we might consider to be the data of the initial-value ODE problem: the function $f(t,u)$, and the initial value $u_0$. It's easier to discuss perturbations to numbers than to functions, so we will focus on the effect of $u_0$ on the solution, using the following theorem that we give without proof. Happily, its conditions are identical to those in the [existence--uniqueness theorem](theorem-existunique).

(theorem-depIC)=

````{prf:theorem} IC dependence
If the derivative $\frac{\partial f}{\partial u}$ exists and $\left|\frac{\partial f}{\partial u}\right|$ is bounded by a constant $L$ for all $a\le t \le b$ and all $u$, then the solution $u(t;u_0+\delta)$ of $u'=f(t,u)$ with initial condition $u(0)=u_0+\delta$ satisfies
  
```{math}
:label: depIC
\left\|u(t;u_0+\delta)-u(t;u_0)\right\|_\infty \le |\delta| e^{L(b-a)},
```

for all sufficiently small $|\delta|$.
````

```{prf:example} Julia demo
:class: demo
:label: demos-basics-cond
{doc}`demos/basics-cond`
```

Numerical solutions of IVPs have errors, and those errors can be seen as perturbations to the solution. The [theorem](theorem-depIC) gives an upper bound of $e^{L(b-a)}$ on the infinity norm (i.e., pointwise) absolute condition number of the solution with respect to perturbations at an initial time. However, the upper bound may be a terrible overestimate of the actual sensitivity for a particular problem.

In general, solutions can diverge from, converge to, or oscillate around the original trajectory in response to perturbations. We won't fully consider these behaviors and their implications for numerical methods again until a later chapter.

## Exercises

1. ✍ For each IVP, determine whether the problem satisfies the conditions of [the IC dependence theorem](theorem-depIC). If so, determine the smallest possible value for $L$.

    **(a)** $f(t,u) = 3 u,\; 0 \le t \le 1$%,\  -\infty < u < \infty.$

    **(b)** $f(t,u) = -t \sin(u),\; 0 \le t \le 5$%,\ -\infty < u < \infty.$

    **(c)** $f(t,u) = -(1+t^2) u^2,\; 1 \le t \le 3$%,\ -\infty < u < \infty.$

    **(d)** $f(t,u) = \sqrt{u},\; 0 \le t \le 1$
  
    ````{only} solutions
    ````

2. ⌨ Solve each IVP in the preceding problem with `DifferentialEquations.solve`, and make a plot of the solution.

    ````{only} solutions
    ````

3. ✍ Use an integrating factor to find the solution of each problem in analytic form.

    **(a)** $u' = -t u,\ 0 \le t \le 5,\ u(0) = 2.$

    **(b)** $u' - 3 u = e^{-2t},\ 0 \le t \le 1,\  u(0) = 5.$

    ````{only} solutions
    ````

4. ✍ Consider the IVP $u'=u^2$, $u(0)=\alpha$.

    **(a)** Does [the existence--uniqueness theorem](theorem-existunique) apply to this problem?

    **(b)** Show that $u(t) = \alpha/(1-\alpha t)$ is a solution of the IVP.

    **(c)** Does this solution necessarily exist for all $t\in[0,1]$?
  
    ````{only} solutions
    ````

    ```{index} logistic equation
    ```

5. ⌨ Using `solve`, compute solutions $x(t)$ to the logistic equation with harvesting,

    ```{math}
    x' = k (S-x)(x-M), \qquad 0\le t \le 10,
    ```

    with $k=S=1$ and $M=0.25$, for the initial conditions $x(0)=0.9M$, $1.1M$, $1.5M$, $0.9S$, $1.1S$, $3S$. Show all the solutions together on one plot. (Note: One of the solutions will throw a warning and fail to reach $t=10$.

    %You can plot it anyway, and make it look better with {term}`ylim`; use \verb+ylim([0,3.5])+ at the end to create a reasonable scale for the vertical axis.)

    ````{only} solutions
    k = 1; S = 1; M = .25;

    f = @(t,x) k*(S-x)*(x-M)
    ts = [0 10];

    for x0 = [.9*M 1.1*M 1.5*M 0.9*S 1.1*S 3*S]
        [t,x] = ode45(f,ts,x0);
        plot(t,x)
        hold on
    end
    ylim([0 3.5])
    ````

6. ⌨ Using `solve`, solve the IVP $u'=u\cos(u) + \cos(4t)$, $0\le t \le 10$, $u(0)=u_0$, for $u_0 = -2,-1.5,-1,\ldots,1.5,2$. Plot all the solutions on a single graph. You should find that they all settle into one of two periodic oscillations. To two digits of accuracy, find the value of $u_0$ in $(-2,2)$ at which the attracting solution changes.

    ````{only} solutions
    ````

    (problem-caffeine)=
7. ⌨ Experimental evidence (see {cite}`newtonPlasmaSalivary1981`) shows that a 300 mg oral dose of caffeine, such as might be found in a large mug of drip-brewed coffee, creates a concentration of about 8 $\mu{\rm g}$/mL in blood plasma. This boost is followed by first-order kinetics with a half-life of about 6 hours (although this rate can vary a great deal from person to person). We can model the caffeine concentration due to one drink taken over half an hour via

    ```{math}
      x'(t) = -kx + C(t),\quad x(0)=0,
    ```

    where $k=\log(2)/6$ and

    ```{math}
      C(t) =
      \begin{cases}
        16, & 0\le t \le 0.5, \\
        0, & t > 0.5.
      \end{cases}
    ```

    Use `solve` to make a plot of the caffeine concentration for 12 hours. Then change $k=\log(2)/8$ (half-life of 8 hours) and plot the solution again.

    ````{only} solutions
    %% (a)
    k = log(2)/6;
    C = @(t) 16*double(t<0.5);
    f = @(t,x) -k*x + C(t);
    [t,x] = ode45(f,[0 12],0);
    plot(t,x)

    %% (b)
    k = log(2)/8;
    C = @(t) 16*double(t<0.5);
    f = @(t,x) -k*x + C(t);
    [t,x] = ode45(f,[0 12],0);
    plot(t,x)
    ````

8. ⌨ A reasonable model of the velocity $v(t)$ of a skydiver is

   $$
   \frac{dv}{dt} = -g + \frac{k}{m}v^2,  \qquad v(0)=0,
   $$

   where $g=9.8 \text{ m/sec}^2$ is gravitational acceleration, $m$ is the mass of the skydiver with parachute, and $k$ quantifies the effect of air resistance. At the US Air Force Academy, a training jump starts at about 1200 m and has $k=0.4875$ for $t<13$ and $k=29.16$ or $t\ge 13$. (This is an oversimplification; see~{cite}`meadeDifferentialEquations1999`.) Find the time at which the skydiver reaches the ground. Keep in mind that the distance fallen up to time $t$ is  $\displaystyle\int_0^t v(s)\, ds$, and use the output form of `solve` as shown in {doc}`demos/basics-usage`.

    ````{only} solutions
    ````
