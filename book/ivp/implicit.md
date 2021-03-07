# Implementation of multistep methods

```{index} multistep method; implementation of
```

We now consider some of the practical issues that arise when multistep formulas are used to solve IVPs. Implementation of the explicit case is relatively straightforward. In what follows we use boldface for the vector form of the problem, $\mathbf{u}'=\mathbf{f}(t,\mathbf{u})$. For instance, the explicit AB4 method is defined by the formula

```{math}
:label: ab4
\mathbf{u}_{i+1} = \mathbf{u}_i + \frac{h}{24} ( 55\mathbf{f}_i - 59 \mathbf{f}_{i-1} +37\mathbf{f}_{i-2} -9\mathbf{f}_{i-3}), \quad i=3,\ldots,n-1.
```

```{prf:example} Julia demo
:class: demo
:label: demos-implicit-ab4
{doc}`demos/implicit-ab4`
```

{numref}`Function {number}<function-ab4>` shows a basic implementation of this formula. Observe that {numref}`Function {number}<function-rk4>` is used to find the starting values $\mathbf{u}_1,\mathbf{u}_2,\mathbf{u}_3$ that are needed before the iteration formula takes over. As far as RK4 is concerned, it needs to solve the IVP over the time interval $a \le t \le a+3h$, using a step size $h$ (the same step size as in the AB4 iteration). These values are then used by {numref}`Function {number}<function-ab4>` to find
$\mathbf{f}_0,\ldots,\mathbf{f}_3$ and get the main iteration started.

For each value of $i$ the formula uses the four most recently known values of the solution's derivative in order to advance by one step. In {numref}`Function {number}<function-ab4>` only these values of $\mathbf{f}$ are stored, and a matrix-vector product is used for the linear combination implied in {eq}`ab4`:

```{math}
  :label: ab4mv
  \mathbf{u}_{i+1} = \mathbf{u}_i + h
  \begin{bmatrix}
    \mathbf{f}_i & \mathbf{f}_{i-1} & \mathbf{f}_{i-2} & \mathbf{f}_{i-3}
  \end{bmatrix}
  \begin{bmatrix}
    55/24 \\[1mm] -59/24 \\[1mm] 37/24 \\[1mm]-9/24
  \end{bmatrix}.
```

```{index} generating polynomials
```

We have distributed the factor of $1/24$ in order to point out that the $4\times 1$ constant vector is just the vector of coefficients of the generating polynomial $\sigma(z)$ from {eq}`sigma`. At the start of an iteration, the value of $\mathbf{f}$ at the most recent solution step is unknown, so a call is made to evaluate it, and the other columns are shifted to the right (i.e., into the past).

(function-ab4)=

````{proof:function} ab4
**4th-order Adams--Bashforth formula for an IVP.**

```{code-block} julia
:lineno-start: 1
"""
ab4(ivp,n)

Apply the Adams-Bashforth 4th order method to solve the given IVP
using `n` time steps.
"""
function ab4(ivp,n)
    # Time discretization.
    a,b = ivp.tspan
    h = (b-a)/n
    t = [ a + i*h for i in 0:n ]

    # Constants in the AB4 method.
    k = 4;    sigma = [55, -59, 37, -9]/24;

    # Find starting values by RK4.
    u = fill(float(ivp.u0),n+1)
    rkivp = ODEProblem(ivp.f,ivp.u0,(a,a+(k-1)*h),ivp.p)
    ts,us = rk4(rkivp,k-1)
    u[1:k] = us[1:k]

    # Compute history of u' values, from newest to oldest.
    f = [ ivp.f(u[k-i],ivp.p,t[k-i]) for i in 1:k-1  ]

    # Time stepping.
    for i in k:n
      f = [ ivp.f(u[i],ivp.p,t[i]), f[1:k-1]... ]   # new value of du/dt
      u[i+1] = u[i] + h*sum(f[j]*sigma[j] for j in 1:k)  # advance a step
    end
    return t,u
end
```
````

## Implicit methods

(function-am2)=

````{proof:function} am2
**2nd-order Adams--Moulton (trapezoid) formula for an IVP.**

```{code-block} julia
:lineno-start: 1
"""
am2(ivp,n)

Apply the Adams-Moulton 2nd order method to solve given IVP using
`n` time steps.
"""
function am2(ivp,n)
    # Time discretization.
    a,b = ivp.tspan
    h = (b-a)/n
    t = [ a + i*h for i in 0:n ]

     # Initialize output.
     u = fill(float(ivp.u0),n+1)

    # Time stepping.
    for i in 1:n
        # Data that does not depend on the new value.
        known = u[i] + h/2*ivp.f(u[i],ivp.p,t[i])
        # Find a root for the new value.
        F = z -> z .- h/2*ivp.f(z,ivp.p,t[i+1]) .- known
        unew = levenberg(F,known)
        u[i+1] = unew[end]
    end
    return t,u
end
```
````

The implementation of an implicit multistep method is a bit more involved. Consider the second-order implicit formula AM2, also known as the trapezoid method. To advance from step $i$ to $i+1$, we need to solve

```{math}
  :label: AM2solve
  \mathbf{z} - \mathbf{u}_i - \tfrac{1}{2} h \bigl[ \mathbf{f}(t_i,\mathbf{u}_i) + f(t_{i+1},\mathbf{z}) \bigr] = 0
```

for $\mathbf{z}$, and then set $\mathbf{u}_{i+1}=\mathbf{z}$. This equation takes the form $\mathbf{g}(\mathbf{z})=\boldsymbol{0}$, so we have a rootfinding problem as in {doc}`../nonlineqn/overview`. An implementation of AM2 using {numref}`Function {number}<function-levenberg>` from {doc}`../nonlineqn/quasinewton` is shown in {numref}`Function {number}<function-am2>`. 

It defines a nested function called `trapzero` that evaluates the left-hand side of {eq}`AM2solve`, given any value of $\mathbf{z}$. The time stepping iteration calls {numref}`Function {number}<function-levenberg>` at each step, starting from the value $\mathbf{u}_i+\tfrac{1}{2}h\mathbf{f}_i$ that is halfway between $\mathbf{u}_i$ and the Euler step $\mathbf{u}_i+h\mathbf{f}_i$. A robust code would have to intercept the case where {numref}`Function {number}<function-levenberg>` fails to converge, but we have ignored this issue for the sake of simplicity.

## Stiff problems

```{prf:example} Julia demo
:class: demo
:label: demos-implicit-stiff
{doc}`demos/implicit-stiff`
```

At each time step in {numref}`Function {number}<function-am2>` (or any implicit IVP solver), a rootfinding iteration of unknown length is needed. This fact makes the cost of an implicit method much greater on a per-step basis than for an explicit one. Given this drawback, you are justified to wonder whether implicit methods are ever competitive! The answer is emphatically yes, as {prf:ref}`demos-implicit-stiff` demonstrates.

```{index} stiff differential equation
```

Although the result of {prf:ref}`demos-implicit-stiff` may seem strange, there is no contradiction: a fourth-order explicit formula is indeed more accurate than a second-order implicit one, in the limit $h\to 0$. But there is another limit to consider, $t\to \infty$ with $h$ fixed, and in this one the implicit method wins. Such problems are called {term}`stiff`. A complete mathematical description will wait for a later chapter, but a sure sign of stiffness is the presence of phenomena on widely different time scales. In the example, there is "slow time," where the solution changes very little, and "fast time," when it suddenly jumps from zero to one. For stiff problems, implicit methods are usually preferred, because they can take far fewer steps than an explicit method, more than offsetting the extra work required per step.

## Adaptivity

As with RK methods, we can run two time stepping methods simultaneously in order to estimate the error and adjust the step size accordingly. For example, we could pair AB3 with AB4 as practically no cost, because the methods differ only in how they include known information from the recent past. The more accurate AB4 value should allow an accurate estimate of the local error in the AB3 value, and so on.

Because multistep methods rely on the solution history, though, changing the step size is more complicated than for RK methods. If $h$ is changed, then the historical values $\mathbf{u}_{i-1},\mathbf{u}_{i-2}\ldots$ and $\mathbf{f}_{i-1},\mathbf{f}_{i-2}\ldots$ are no longer given at the right moments in time to apply the iteration formula. A typical remedy is to use interpolation to re-evaluate the historical values at the appropriate times. The details are important but not especially illuminating, and we do not give them here.

## Exercises

% Must stay as #1
(problem-ab4tests)=

1. ⌨ For each IVP, use {numref}`Function {number}<function-ab4>` to find the solution over the indicated time interval for $n=250$. Plot the computed solution $(t_i,u_i)$ for $i=0,\ldots,n$, and separately plot the error $\bigl(t_i,u_i-\hat{u}(t_i)\bigr)$.

    **(a)** $u' = -2t u, \ 0 \le t \le 2, \ u(0) = 2;\  \hat{u}(t) = 2e^{-t^2}$

    **(b)** $u' = u + t, \ 0 \le t \le 1, \ u(0) = 2;\  \hat{u}(t) = 1-t+e^t$

    **(c)** $u' = x^2/[u(1+x^3)],\ 0 \le x \le 3, \ u(0) =1;\ \hat{u}(x) =[1+(2/3)\ln (1+x^3)]^{1/2}$

    **(d)** $u''+ 9u = 9t, \: 0< t< 2\pi, \: u(0) =1,\: u'(0) = 1; \: \hat{u}(t) = t+\cos (3t)$

    **(e)** $u''+ 9u = \sin(2t), \: 0< t< 2\pi, \: u(0) =2,\: u'(0) = 1$;
    $\quad \hat{u}(t) = (1/5) \sin(3t) + 2 \cos (3t)+ (1/5) \sin (2t)$

    **(f)** $u''- 9u = 9t \: 0< t< 1, \: u(0) =2,\: u'(0) = -1; \: \hat{u}(t) = e^{3t} + e^{-3t}-t$

    **(g)** $u''+ 4u'+ 4u = t, \: 0< t< 4, \: u(0) =1,\: u'(0) = 3/4; \: \hat{u}(t) = (3t+5/4)e^{-2t} + (t-1)/4$

    **(h)** $x^2 u'' +5xu' + 4u = 0,\: 1<x<e^2, \: u(1) =1, \: u'(1) = -1; \: \hat{u}(x) = x^{-2}( 1 + \ln x)$

    **(i)** $2 x^2 u'' +3xu' - u = 0,\: 1<x<16, \: u(1) =4, \: u'(1) = -1$;
    $\quad \hat{u}(x) = 2(x^{1/2} + x^{-1})$

    **(j)** $x^2 u'' -xu' + 2u = 0,\: 1<x<e^{\pi}, \: u(1) =3, \: u'(1) = 4$;
    $\quad \hat{u}(x) = x \left[ 3 \cos \left( \ln x \right)+\sin \left( \ln x \right) \right]$

    ````{only} solutions
    ````

    % Must stay as #2
    (problem-ab4converge)=

2. ⌨ For each IVP in the preceding problem, use {numref}`Function {number}<function-ab4>` for $n=10\cdot2^d$ and $d=1,\ldots,10$. Make a log-log convergence plot for the final time error $|u_n-\hat{u}(t_n)|$ versus $n$, and add a straight line indicating fourth-order convergence.

    ````{only} solutions
    ````

3. ⌨  Line 35 of {numref}`Function {number}<function-ab4>` reads

    ``` julia
    u(:,i+1) = u(:,i) + h*(f*sigma);
    ```

    Explain carefully why this is preferable to

    ``` julia
    u(:,i+1) = u(:,i) + h*f*sigma;
    ```

4. ⌨ Repeat [exercise 1](problem-ab4tests) above  using {numref}`Function {number}<function-am2>`.

    ````{only} solutions
    ````

5. ⌨  Repeat [exercise 2 above](problem-ab4converge)  using {numref}`Function {number}<function-am2>` and comparing to second-order rather than fourth-order convergence.

    ````{only} solutions
    ````

6. ⌨ Using {numref}`Function {number}<function-am2>` as a model, write a function `bd2` that applies the BD2 method to solve an IVP. Test the convergence of your function on one of the IVPs in [exercise 1](problem-ab4tests) above.

    ````{only} solutions
    ````

7. ⌨ For numerical purposes, the exact solution of the IVP in {prf:ref}`demos-implicit-stiff` satisfies $\hat{u}(400)=1$.

    **(a)** Use {numref}`Function {number}<function-ab4>` with $n=200,400,600,\ldots,2000$ and make a log-log convergence plot of the error $|u_n-1|$ as a function of $n$.

    **(b)** Repeat part~(a) using {numref}`Function {number}<function-am2>`.
  
    ````{only} solutions
    %% (a)
    n_ = 200:200:2000;  err_ = [];
    for n = n_
        [t,u] = ab4(@(t,u) u^2-u^3,[0 400],0.005,n);
        err_ = [err_;abs(u(end)-1)];
    end
    clf, loglog(n_,err_,'-o')
    %% (b)
    err_ = [];
    for n = n_
        [t,u] = am2(@(t,u) u^2-u^3,[0 400],0.005,n);
        err_ = [err_;abs(u(end)-1)];
    end
    hold on, loglog(n_,err_,'-o')
    ````

    (problem-ivpimag)=

8. Consider the IVP
  
    ```{math}
    \mathbf{u}'(t) = \mathbf{A} \mathbf{u}(t), \quad \mathbf{A}=
    \begin{bmatrix}
      0&-4\\4&0
    \end{bmatrix}, \quad \mathbf{u}(0) =
    \begin{bmatrix}
      1\\0
    \end{bmatrix}.
    ```

    **(a)** ✍ Define $E(t) = \bigl\|\mathbf{u}(t)\bigr\|_2^2$. Show that $E(t)$ is constant. (Differentiate $u^Tu$ with respect to time and show that it simplifies to zero.)

    **(b)** ⌨ Use {numref}`Function {number}<function-ab4>` to solve the IVP for $t\in[0,20]$ with $n=100$ and $n=150$. On a single graph using a log scale on the $y$-axis, plot $|E(t)-E(0)|$ versus time for both solutions. You should see exponential growth in time.

    **(c)** ⌨ Repeat part~(b) with $n=400$ and $n=600$, but use a linear scale on the $y$-axis. Now you should see only linear growth of $|E(t)-E(0)|$.
  
    ````{only} solutions
    %%
    %% (a)
    % $$dE/dt = 2u^T\dot{u} = 2u^T(Au) = 2(-4u_1u_2+4u_2u_1)=0.$$
    %% (b)
    A = [0 -4; 4 0];
    clf
    for n = [100 150 200]
        [t,u] = ab4(@(t,u) A*u,[0 20],[1;0],n);
        E = sum(u.^2,2);
        semilogy(t,abs(E-E(1))), hold on
    end
    %%
    A = [0 -4; 4 0];
    clf
    for n = [200 400 800]
        [t,u] = ab4(@(t,u) A*u,[0 20],[1;0],n);
        E = sum(u.^2,2);
        plot(t,E-E(1)), hold on
    end
    ````

9. ⌨ **(a)** Modify {numref}`Function {number}<function-ab4>` to implement the AB2 method.

    **(b)** Repeat part (b) of the preceding exercise, using AB2 in place of AB4.

    **(c)** Repeat part (c) of the preceding exercise, using AB2 in place of AB4.
  

    ````{only} solutions
    %% (a)
    function [t,u] = ab2(dudt,tspan,u0,n)

    % Discretize time.
    a = tspan(1);  b = tspan(2);
    h = (b-a)/n;
    t = tspan(1) + (0:n)'*h;

    k = 2;
    sigma = [3; -1]/2;

    u = zeros(length(u0),n+1);
    [ts,us] = rk4(dudt,[a, a+(k-1)*h],u0,k-1);
    u(:,1:k) = us(1:k,:)';

    f = zeros(length(u0),k);
    for i = 1:k-1
      f(:,k-i) = dudt(t(i),u(:,i));
    end

    % Time stepping.
    for i = k:n
      f = [dudt(t(i),u(:,i)), f(:,1:k-1)];   % new value of du/dt
      u(:,i+1) = u(:,i) + h*(f*sigma);       % advance one step
    end

    u = u.';   % conform to MATLAB output convention
    %% (b)
    A = [0 -4; 4 0];
    clf
    for n = [100 150 200]
        [t,u] = ab2(@(t,u) A*u,[0 20],[1;0],n);
        E = sum(u.^2,2);
        semilogy(t,abs(E-E(1))), hold on
    end
    %%
    A = [0 -4; 4 0];
    clf
    for n = [200 400 800]
        [t,u] = ab2(@(t,u) A*u,[0 20],[1;0],n);
        E = sum(u.^2,2);
        plot(t,E-E(1)), hold on
    end
    ````
