# Adaptive integration

```{index} numerical integration
```

````{prf:example} Julia demo
:class: demo
:label: demos-adapt-motive
{doc}`demos/adapt-motive`
````

To this point, we have used only equally spaced nodes to compute integrals. Yet there are problems in which non-uniformly distributed nodes would clearly be more appropriate, as demonstrated in {prf:ref}`demos-adapt-motive`.

We would like an algorithm that automatically detects and reacts to a situation like that in {prf:ref}`demos-adapt-motive`, a trait known as **adaptivity**.

## Error estimation

```{index} Simpson formula
```

Ideally, we would like to make adaptation decisions based on the error of the integration. Knowing the error exactly would be equivalent to knowing the exact answer, but we can estimate it using the extrapolation technique of {doc}`integration`. Consider the Simpson formula {eq}`extraplevel1` resulting from one level of extrapolation from trapezoid estimates:

```{math}
  :label: extraplevel1repeat
  S_f(2n) = \frac{1}{3} \Bigl[ 4 T_f(2n) - T_f(n) \Bigr].
```

We expect this method to be fourth-order accurate, i.e.,

```{math}
  \int_a^b f(x)\, dx = S_f(2n) + O(n^{-4}),
```

We can further extrapolate to sixth-order accuracy by {eq}`nc-sixth`,

```{math}
  :label: extraplevel2repeat
  R_f(4n) = \frac{1}{15} \Bigl[ 16 S_f(4n) - S_f(2n) \Bigr].
```

By virtue of higher order of accuracy, $R_f(4n)$ should be more accurate than $S_f(4n)$. Hence a decent estimate of the error in the better of the two Simpson values is

```{math}
  :label: adapterr
  E = R_f(4n) - S_f(4n) = \frac{S_f(4n) - S_f(2n)}{15}.
```

## Divide and conquer

If $|E|$ is judged to be acceptably small, we are done. This judgment takes some care. For instance, suppose the exact integral is $10^{20}$.  Requiring $|E| < \delta\ll 1$ would be fruitless in double precision, since it would require more than 20 accurate digits. Hence checking the absolute size of the error alone is not appropriate. Conversely, consider the integral

```{math}
  \int_{10^{-6}}^{2\pi} 2 \sin x\, dx \approx -10^{-12}.
```

We are likely to sample values of the integrand that are larger than, say, $1/2$ in absolute value, so obtaining this very small result has to rely on subtractive cancellation. We cannot hope for more than 4-5 accurate digits, so a strict test of the relative error is also not recommended.\footnote{In other words, we can have an error that is small relative to the data (the integrand), which is $O(1)$, but not relative to the answer itself.} Typically we use both relative and absolute error, stopping when either one is considered small enough. Algebraically, the test is

```{math}
  :label: absreltolerance
  |E| < \delta_a + \delta_r |S_f(n)|,
```

where $\delta_a$ and $\delta_r$ are given absolute and relative error tolerances, respectively.

When $|E|$ fails to meet {eq}`absreltolerance`, we bisect the interval $[a,b]$ to exploit the identity

```{math}
  \int_a^b f(x)\, dx = \int_a^{(a+b)/2} f(x)\, dx + \int_{(a+b)/2}^b f(x)\, dx,
```

and independently compute estimates to each of the half-length integrals. Each of these half-sized computations recursively applies Simpson's formula and the error estimation criterion, making further bisections as necessary. Such an approach is called **divide and conquer** in computer science: recursively split the problem into easier pieces and glue the results together.

## Implementation

It is typical to use just the minimal formula $S_f(4)$ and its error estimate $E$ to make decisions about adaptivity. A computation of $S_f(4)$ requires three trapezoid estimates $T_f(1)$, $T_f(2)$, and $T_f(4)$. As observed in {eq}`nc-doubling` and {prf:ref}`demos-int-extrap`, the five integrand evaluations in $T_f(4)$ are sufficient to compute all of these values. There is one further exploitation of node locations to be found. For simplicity, assume $[a,b]=[0,1]$. The five nodes used in $T_f(4)$ are

```{math}
0, \quad \frac{1}{4}, \quad  \frac{1}{2}, \quad  \frac{3}{4}, \quad 1.
```

If we bisect the interval and compute $T_f(4)$ on the subinterval $[0,1/2]$, we use the nodes

```{math}
0, \quad \frac{1}{8}, \quad  \frac{1}{4}, \quad  \frac{3}{8}, \quad \frac{1}{2}.
```

Only the second and fourth nodes are new. The same is true on the subinterval $[1/2,1]$, and for every recursive bisection.

````{prf:example} Julia demo
:class: demo
:label: demos-adapt-usage
{doc}`demos/adapt-usage`
````

{numref}`Function {number}<function-intadapt>` shows how to exploit this structure. The nested function `do_integral` does all of the work. It expects to receive the three nodes and integrand values that it shares with the level above. It adds the two new nodes and uses the set of all five to compute three trapezoid estimates with $n=1$, $n=2$, and $n=4$, using the updating formula {eq}`nc-doubling` twice. It goes on to find the two Simpson approximations and to estimate the error in the better one by {eq}`adapterr`.

If the error estimate passes the test {eq}`absreltolerance`, the better Simpson value is returned as the integral over the given interval. Otherwise, the interval is bisected, the two pieces computed using recursive calls, and those results are added to give the complete integral.

(function-intadapt)=

````{proof:function} intadapt
**Adaptive integration with error estimation**

```{code-block} julia
:lineno-start: 1
"""
intadapt(f,a,b,tol)

Do adaptive integration to estimate the integral of `f` over
[`a`,`b`] to desired error tolerance `tol`. Returns estimate and a
vector of evaluation nodes used.
"""
function intadapt(f,a,b,tol)
    # Use error estimation and recursive bisection.
    function do_integral(a,fa,b,fb,m,fm,tol)
        # These are the two new nodes and their f-values.
        xl = (a+m)/2;  fl = f(xl);
        xr = (m+b)/2;  fr = f(xr);
        t = [a,xl,m,xr,b]          # all 5 nodes at this level

        # Compute the trapezoid values iteratively.
        h = (b-a)
        T = [0.,0.,0.]
        T[1] = h*(fa+fb)/2
        T[2] = T[1]/2 + (h/2)*fm
        T[3] = T[2]/2 + (h/4)*(fl+fr)

        S = (4*T[2:3]-T[1:2]) / 3      # Simpson values
        E = (S[2]-S[1]) / 15           # error estimate

        if abs(E) < tol*(1+abs(S[2]))  # acceptable error?
            Q = S[2]                   # yes--done
        else
            # Error is too large--bisect and recurse.
            QL,tL = do_integral(a,fa,m,fm,xl,fl,tol)
            QR,tR = do_integral(m,fm,b,fb,xr,fr,tol)
            Q = QL + QR
            t = [tL;tR[2:end]]   # merge the nodes w/o duplicate
        end
        return Q,t
    end

    m = (b+a)/2
    Q,t = do_integral(a,f(a),b,f(b),m,f(m),tol)
    return Q,t
end
```
````

Although adaptivity and the error estimation that goes with it can be very powerful, they come at some cost. The error estimation cannot be universally perfect, so sometimes the answer will not be as accurate as requested (underestimation) and sometimes the function will be evaluated more times than necessary (overestimation). Subtle problems may arise when the integral is a step within a larger computation (see [this exercise below](problem-adapt-int-nonsmooth)).

## Exercises

(problem-adaptquadtests)=
% must be kept as #1

1. ⌨ For each integral below, use {numref}`Function {number}<function-intadapt>` with error tolerance $10^{-2}$, $10^{-3}$, \ldots, $10^{-12}$. Make a table of errors and the number of integrand evaluation nodes used, and use a convergence plot as in {prf:ref}`demos-adapt-usage` to compare to fourth-order accuracy. (These integrals were taken from {cite}`baileyComparisonThree2005`.)

    **(a)** $\displaystyle \int_0^1 x\log(1+x)\, dx = \frac{1}{4}$

    **(b)** $\displaystyle \int_0^1 x^2 \tan^{-1}x\, dx = \frac{\pi-2+2\log 2}{12}$

    **(c)** $\displaystyle \int_0^{\pi/2}e^x \cos x\, dx = \frac{e^{\pi/2}-1}{2}$

    **(d)** $\displaystyle \int_{0}^1 \sqrt{x} \log(x) \, dx = -\frac{4}{9}$ (Note: Although the integrand has the limiting value zero as $x\to 0$, you have to implement the function carefully to return zero as the value of $f(0)$, or start the integral at $x=\macheps$.)

    **(d)** $\displaystyle \int_0^1 \sqrt{1-x^2}\, dx = \frac{\pi}{4}$

    ````{only} solutions
    f = {
        @(x) x.*log(1+x);
        @(x) x.^2.*atan(x);
        @(x) exp(x).*cos(x);
        @(x) sqrt(x).*log(x);
        @(x) sqrt(1-x.^2);
        };

    a = [ 0;0;0;eps;0];
    b = [1;1;pi/2;1;1];
    I = [.25;(pi-2+2*log(2))/12;(exp(pi/2)-1)/2;-4/9;pi/4];

    tol_ = 10.^(-1:-1:-14)';
    err_ = [];
    num_ = [];

    for tol = tol_'
        err = [];  num = [];
        for i = 1:length(f)
            [Q,t] = intadapt(f{i},a(i),b(i),tol);
            err = [err, I(i) - Q];
            num = [num, length(t)];
        end
        err_ = [err_; err];
        num_ = [num_; num];
    end

    loglog(num_,abs(err_),'o-')
    ````

2. ⌨ For each integral below: (i) use `quadgk` to find the value to at least 12 digits; (ii) use {numref}`Function {number}<function-intadapt>` to evaluate the integral to a tolerance of $10^{-8}$; (iii) compute the absolute error and the number of nodes used; (iv) use the $O(h^2)$ term in the Euler--Maclaurin formula {eq}`eulermaclaurin` to estimate how many nodes are required by the fixed-stepsize trapezoidal formula to reach an absolute error of $10^{-8}$.

    **(a)** $\displaystyle \int_{0.1}^3 \operatorname{sech}(\sin(1/x))\, d x$

    **(b)** $\rule[2em]{0pt}{0pt} \displaystyle\int_{-0.9}^9 \ln((x+1)^3))\, d x$

    **(c)** $\rule[2em]{0pt}{0pt} \displaystyle\int_{-\pi}^\pi \cos(x^3)\, d x$
  
    ````{only} solutions
    format long
    %% (a)
    % $$ \int_{0.1}^3 {\rm sech}(\sin(1/x))\,dx$$
    %
    %%
    f = @(x) sech( sin(1./x) );
    dfdx = @(x) +sech(sin(1./x)).*tanh(sin(1./x)).*cos(1./x)./x.^2;
    [Q,x] = intadapt(f,0.1,3,1e-8);
    Q
    n_adapt = length(x)-1
    I = integral(f,0.1,3,'abstol',1e-13,'reltol',1e-13);
    abs_err = abs(Q-I)
    rel_err = abs_err/abs(I)

    %%
    h_trap = sqrt( 12e-8 / abs(dfdx(1)-dfdx(0.1)) )
    n_trap = ceil( (3-0.1)/h_trap )


    %% (b)
    % $$\int_{-0.9}^9 \ln((x+1)^3))\, d x$$
    %%
    f = @(x) log( (x+1).^3 );
    dfdx = @(x) 3./(x+1);
    [Q,x] = intadapt(f,-0.9,9,1e-8);
    Q
    n_adapt = length(x)-1

    %%
    h_trap = sqrt( 12e-8 / abs(dfdx(9)-dfdx(-0.9)) )
    n_trap = ceil( (9+0.9)/h_trap )


    %% (c)
    % $$\int_{-\pi}^\pi \cos(x^3)\,d x$$
    %%
    f = @(x) cos( x.^3 );
    dfdx = @(x) -sin(x.^3).*3*x.^2;
    [Q,x] = intadapt(f,-pi,pi,1e-8);
    Q
    n_adapt = length(x)-1

    %%
    h_trap = sqrt( 12e-8 / abs(dfdx(pi)-dfdx(-pi)) )
    n_trap = ceil( (pi+pi)/h_trap )
    ````

3. ⌨ An integral such as $\displaystyle \int_0^1 x^{-\gamma}\, dx$ for $\gamma>0$, in which the integrand blows up at one or both ends, is known as an **improper** integral. It has a finite value if $\gamma<1$, despite the singularity. One way to deal with the problem of the infinite value for $f(t_0)$ is to replace the lower limit with a small number $\epsilon$. Using {numref}`Function {number}<function-intadapt>` with a small tolerance, make a log--log plot of the error as a function of $\epsilon$ for $\epsilon=10^{-15},10^{-16},\ldots,10^{-45}$. (A more robust way to handle improper integrals is discussed in a later chapter.)

    ````{only} solutions
    ````

4. ⌨ A curious consequence of our logic in {numref}`Function {number}<function-intadapt>` is that the algorithm uses what we believe to be a more accurate, sixth-order answer only for estimating error; the returned value is the supposedly less accurate $S_f(2n)$. The practice of returning the extrapolated $R_f(4n)$ instead is called {index}`local extrapolation` *local extrapolation*. Modify {numref}`Function {number}<function-intadapt>` to use local extrapolation and repeat problem 1 above. Is the convergence more like 4th order or 6th order?

    ````{only} solutions
    ````

5. ⌨ The **sine integral function** is defined by
  
    ```{math}
    \operatorname{Si}(x) = \int_0^x \frac{\sin z}{z}\, dz.
    ```

    Use {numref}`Function {number}<function-intadapt>` to plot Si over the interval $[1,10]$. Note: You will need to replace the lower bound of integration by $\macheps$.

    ````{only} solutions
    ````

    (problem-adapt-int-nonsmooth)=

6. ⌨  Adaptive integration can have subtle drawbacks. This exercise is based on the **error function**, a smooth function defined as
  
    ```{math}
    \operatorname{erf}(x) = \frac{2}{\pi}\int_0^x e^{-s^2}\,ds.
    ```

    **(a)** Define a function $g$ that approximates erf by applying {numref}`Function {number}<function-trapezoid>` with $n=100$. Make a plot of the error $g(x)-\operatorname{erf}(x)$ at 300 points in the interval $[0,3]$.

    **(b)** Define another approximation $h$ that applies {numref}`Function {number}<function-intadapt>` with error tolerance $10^{-7}$. Plot the error in $h$ as in part~(a). Why does it look so different from the previous case?

    **(c)** Suppose you wished to find $x$ such that $\operatorname{erf}(x) = .95$ by using rootfinding on one of your two approximations. Which would be preferable?
