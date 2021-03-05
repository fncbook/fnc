# Adaptive Runge--Kutta

```{index} Runge--Kutta method
```

The derivation and analysis of methods for initial-value problems usually assumes a fixed step size $h$. While the error behavior $O(h^p)$ is guaranteed by the [convergence theorem](theorem-onestepGTE) as $h\rightarrow 0$, this bound comes with an unknowable constant, and it is not very useful as a guide to the numerical value of the error at any particular value of $h$. Furthermore, as we saw in {doc}`../localapprox/adaptive` with numerical integration, in many problems a fixed value of $h$ throughout $a\le t \le b$ is far from the most efficient strategy.

In response we will employ the basic strategy of {doc}`../localapprox/adaptive`: adapt the time step size in order to reach an accuracy goal, as measured by an error estimate formed from computing multiple approximations. The details are quite different, however.

## Error estimation

Suppose that, starting from a given value $u_i$ and using a step size $h$, we run one step of *two* RK methods simultaneously: one method with order $p$, producing $u_{i+1}$, and the other method with order $p+1$, producing $\tilde{u}_{i+1}$. In most circumstances, we can expect that $\tilde{\mathbf{u}}_{i+1}$ is a much better approximation to the solution than $\mathbf{u}_{i+1}$ is. So it seems reasonable to use $E_i(h)=|\tilde{\mathbf{u}}_{i+1}-\mathbf{u}_{i+1}|$ (in the vector case, a norm) as an estimate of the actual local error made by the $p$th order method.

If our goal is to keep error less than some predetermined value $\epsilon$, we could decide to accept the new solution value if $E_i<\epsilon$ and otherwise reject it. (Even though the estimate $E_i$ is meant to go with the *less* accurate proposed value $\mathbf{u}_{i+1}$, it's hard to resist the temptation to keep the more accurate value $\tilde{\mathbf{u}}_{i+1}$ instead, and this is common in practice.)

Regardless of whether $E_i<\epsilon$ and we accept the step, we now ask a question: looking back, what step size *should* we have taken to just meet our error target $\epsilon$? Let's speculate that $E_i(h)\approx C h^{p+1}$ for an unknown constant $C$, given the behavior of local truncation error as $h\rightarrow 0$. If we had used a step size $qh$ for some $q>0$, then trivially, $E_i(qh)\approx C q^{p+1}h^{p+1}$ is what we would expect to get. Our best guess for $q$ would be to set $E_i(qh)\approx \epsilon$, or

```{math}
  :label: adaptRKlocal
  q \approx \left(\frac{\epsilon}{E_i}\right)^{1/(p+1)}.
```

Whether or not we accepted the value proposed for $t=t_{i+1}$, we will adjust the step size to $q h$ for the next attempted step.

Given what we know about the connection between local and global errors, we might instead decide that controlling the normalized contribution to *global* error, which is closer to $E_i(qh)/(q h)$, is more reasonable. Then we end up with

```{math}
  :label: adaptRKglobal
  q \le \left(\frac{\epsilon}{E_i}\right)^{1/p}.
```

Experts have different recommendations about whether to use {eq}`adaptRKlocal` or {eq}`adaptRKglobal`. Even though {eq}`adaptRKglobal` appears to be more in keeping with our assumptions about global errors, modern practice seems to favor {eq}`adaptRKlocal`.

## Embedded formulas

We have derived two useful pieces of information: a reasonable estimate of the actual value of the local (or global) error, and a prediction how the step size will affect that error. Together they can be used to adapt step size and keep errors near some target level. But there remains one more important twist to the story.

```{margin}
Embedded RK formulas are a pair of RK methods whose stages share the same internal $f$ evaluations.
```

At first glance, it would seem that to use (for example) any pair of second- and third-order RK methods to get the $\mathbf{u}_{i+1}$ and $\tilde{\mathbf{u}}_{i+1}$ needed for adaptive error control, we need at least $2+3=5$ evaluations of $f(t,y)$ for each attempted time step.  This is more than double the computational work needed by the second-order method without adaptivity. Fortunately, the marginal cost of adaptation can be substantially reduced by using **embedded Runge--Kutta** Embedded RK formulas are a pair of RK methods whose stages share the same internal $f$ evaluations, combining them differently in order to get estimates of two different orders of accuracy.

A good example of an embedded method is the **Bogacki--Shampine** (BS23) formula, given by the table

```{math}
:label: bs23
\begin{array}{r|cccc}
0                  & \rule{0pt}{2.75ex} &                    &                    &                    \\
\frac{1}{2}        & \frac{1}{2}        & \rule{0pt}{2.75ex} &                    &                    \\
\frac{3}{4}        & 0                  & \frac{3}{4}        & \rule{0pt}{2.75ex} &                    \\
1                 & \frac{2}{9}        & \frac{1}{3}        & \frac{4}{9}        & \rule{0pt}{2.75ex} \\[2pt] \hline
\rule{0pt}{2.75ex} & \frac{2}{9}        & \frac{1}{3}        & \frac{4}{9}        & 0                  \\[2pt] \hline
\rule{0pt}{2.75ex} & \frac{7}{24}       & \frac{1}{4}        & \frac{1}{3}        & \frac{1}{8}
\end{array}
```

The top part of the table describes four stages in the usual RK fashion. The last two rows describe how to construct a second-order estimate $\mathbf{u}_{i+1}$ and a third-order estimate $\tilde{\mathbf{u}}_{i+1}$ by taking different combinations of those stages.

## Implementation

(function-rk23)=

````{proof:function} rk23
**Adaptive IVP solver based on embedded RK formulas.**

```{code-block} julia
:lineno-start: 1
"""
rk23(ivp,tol)

Apply an adaptive embedded RK formula pair to solve given IVP with
estimated error `tol`. Returns a vector of times and a vector of
solution values.
"""
function rk23(ivp,tol)
    # Initialize for the first time step.
    a,b = ivp.tspan
    t = [a]
    u = [float(ivp.u0)];   i = 1;
    h = 0.5*tol^(1/3)
    s1 = ivp.f(ivp.u0,ivp.p,a)

    # Time stepping.
    while t[i] < b
        # Detect underflow of the step size.
        if t[i]+h == t[i]
            @warn "Stepsize too small near t=$(t[i])"
            break  # quit time stepping loop
        end

        # New RK stages.
        s2 = ivp.f( u[i]+(h/2)*s1,   ivp.p, t[i]+h/2   )
        s3 = ivp.f( u[i]+(3*h/4)*s2, ivp.p, t[i]+3*h/4 )
        unew2 = u[i] + h*(2*s1 + 3*s2 + 4*s3)/9   # 2rd order solution
        s4 = ivp.f( unew2, ivp.p, t[i]+h )
        err = h*(-5*s1/72 + s2/12 + s3/9 - s4/8)  # 2nd/3rd difference
        E = norm(err,Inf)                         # error estimate
        maxerr = tol*(1 + norm(u[i],Inf))     # relative/absolute blend

        # Accept the proposed step?
        if E < maxerr     # yes
            push!(t,t[i]+h)
            push!(u,unew2)
            i += 1
            s1 = s4       # use FSAL property
        end

        # Adjust step size.
        q = 0.8*(maxerr/E)^(1/3)   # conservative optimal step factor
        q = min(q,4)               # limit stepsize growth
        h = min(q*h,b-t[i])        # don't step past the end
    end
    return t,u
end
```
````

Our implementation of an embedded second/third order (RK23) code is given in {ref}`function-rk23`. It has a few details that are worth explaining.

First, as in {eq}`absreltolerance`, we use a combination of absolute and relative tolerances to judge the acceptability of a solution value. Second, we have a check whether $t_i+h$ equals $t_i$, which looks odd. This check is purely about roundoff error, because $h$ can become so small that it no longer changes the floating point value of $t_i$. When this happens, it's often a sign that the underlying exact solution has a singularity near $t=t_i$. Third, some adjustments are made to the step size prediction factor $q$. We use a smaller value than {eq}`adaptRKlocal`, to be conservative about the many assumptions that were made to derive it. We also prevent a huge jump in the step size for the same reason. And, we make sure that our final step doesn't take us past the requested end of the domain.

```{proof:example} Julia demo
:class: demo
{doc}`demos/adapt-basic`

{doc}`demos/adapt-sing`
```

Finally, there is some careful programming done to avoid redundant evaluations of $f$. As written in {eq}`bs23`, there seem to be four stages needed to find the paired second- and third-order estimates. This is unfortunate, since there are three-stage formulas of order three. But BS23 has a special property called "first same as last" (FSAL). If the proposed step is accepted, the final stage computed in stepping from $t_i$ to $t_{i+1}$ is identical to the *first* stage needed to step from $t_{i+1}$ to $t_{i+2}$, so in that sense one of the stage evaluations comes at no cost. This detail is addressed in our code.

Often the steps chosen adaptively clearly correspond to identifiable features of the solution. However, there are so-called **stiff problems** in which the time steps seem unreasonably small in relation to the observable behavior of the solution. These problems benefit from a particular type of solver and will be taken up in {doc}`implicit`.

## Exercises

1. ⌨ Using {ref}`function-rk23`, solve $y'' +(1+y')^3 y = 0$ over $ 0 \le t \le 4 \pi$ with the indicated initial conditions. Plot $y(t)$ and $y'(t)$ as a function of $t$ and separately plot the solution curve parametrically in the phase plane—that is, the $\bigl(y(t),y'(t)\bigr)$-plane.

    **(a)** $y(0) = 0.1, \quad y'(0) = 0$

    **(b)** $y(0) = 0.5, \quad y'(0) = 0$

    **(c)** $y(0) = 0.75, \quad y'(0) = 0$

    **(d)** $y(0) = 0.95, \quad y'(0) = 0$

    ````{only} solutions
    ````

2. ⌨ Solve the [caffeine exercise](problem-caffeine) using {ref}`function-rk23` with an error tolerance of $10^{-5}$. Plot the solution so that you can see the individual points. What is the smallest time step taken, and at what time does it occur?

    ````{only} solutions
    k = log(2)/6;
    C = @(t) 16*double(t<0.5);
    f = @(t,x) -k*x + C(t);
    [t,x] = rk23(f,[0 12],0,1e-5);
    plot(t,x,'.-')
    [h,i] = min(diff(t))
    t(i)
    %h =
    %   2.8309e-04
    %i =
    %    16
    %ans =
    %    0.5002
    ````

3. ⌨ Solve the [FitzHugh--Nagumo exercise](problem-fitznag) using {ref}`function-rk23`. Let the error tolerance be $10^{-k}$, increasing the integer $k$ until the graph of the solutions no longer changes. (This illustrates that the error tolerance is a request, not a guarantee!)

    ````{only} solutions
    a = 0.1;
    epsilon = 0.008;
    gamma = 1;

    I_ = [0.05527; 0.05683; 0.0568385; 0.05740];

    clf
    for I = I_'
        f = @(t,v) [ -v(1)*(v(1)-1)*(v(1)-a) - v(2) + I; ...
            epsilon*(v(1)-gamma*v(2)) ];
        [t,v] = rk23(f,[0,600],[0.5;0.1],1e-11);
        plot(t,v(:,1)), grid on, hold on
    end
    % I need a tolerance of about 1e-11 to see no more changes.
    ````

4. ✍ Derive equation {eq}`adaptRKglobal` using the stated assumption about controlling global rather than local error.

    ````{only} solutions

    ````

5. ⌨ Solve the problem $u'=u^2-u^3$, $u(0)=0.001$, $0\le t \le 2000$ and make plots as in {doc}`demos/adapt-basic` that show both the solution and the time steps taken. Does the step size selection seem to be entirely explained by the local variability of the solution?

    ````{only} solutions

    ````
