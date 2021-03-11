# Stability of polynomial interpolation

```{index} stability; of polynomial interpolation
```

With  barycentric interpolation available in the form of {numref}`Function {number}<function-polyinterp>`, we can explore polynomial interpolation using a numerically stable algorithm. Any remaining sensitivity to error is due to the interpolation process itself.

::::{prf:example} Julia demo
:class: demo
:label: demos-stability-equispaced
{doc}`demos/stability-equispaced`
::::

## Runge phenomenon

The disappointing loss of convergence in {ref}`example-equiinterp1` is due to the use of equally spaced nodes. We will examine this effect using the error formula {eq}`interperror` as a guide:

$$
  f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \Phi(x), \qquad \Phi(x) =
  \prod_{i=0}^n (x-t_i).
$$

The $\Phi(x)$ term, which we call the **error indicator function**, can be studied as a function of the nodes only.

::::{prf:example} Julia demo
:class: demo
:label: demos-stability-errfun
{doc}`demos/stability-errfun`
::::

Even though $\Phi\to 0$ at every point in the interval, the exponentially growing gap between the ends and the middle of the interval can ruin the convergence of polynomial interpolation for many choices of $f$.

::::{prf:example} Julia demo
:class: demo
:label: demos-stability-runge
{doc}`demos/stability-runge`
::::

```{index} Runge phenomenon
```
The observation of instability in {prf:ref}`demos-stability-runge` is known as the {term}`Runge phenomenon`. The Runge phenomenon is an instability in the abstract mapping from a function to its polynomial interpolant, manifested when the nodes of the interpolant are equally spaced and the degree of the polynomial increases. We reiterate that the phenomenon is rooted in the convergence theory and not a consequence of the algorithm chosen to implement polynomial interpolation.

Significantly, the convergence observed in {prf:ref}`demos-stability-runge` is stable within a middle portion of the interval. By redistributing the interpolation nodes, we will next sacrifice a little of the convergence in the middle portion in order to improve it enough near the ends to rescue the process globally.

## Chebyshev nodes

```{index} Chebyshev points; second kind
```
The observations above suggest that we might find success by having more nodes near the ends of the interval than in the middle. Though we will not give the details, it turns out that there is a precise asymptotic sense in which this must be done to make polynomial interpolation work over the entire interval. One especially important node family that gives stable convergence for polynomial interpolation is the **Chebyshev points of the second kind**} (also known as Chebyshev extreme points) defined by

:::{math}
:label: chebextreme
 t_k = - \cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n.
:::

These are the projections onto the $x$-axis of $n$ equally spaced points on a unit circle. As such they are densely clustered near the ends of $[-1,1]$, and this feature turns out to overcome the Runge phenomenon.

::::{prf:example} Julia demo
:class: demo
:label: demos-stability-errcheb
{doc}`demos/stability-errcheb`
::::

::::{prf:example} Julia demo
:class: demo
:label: demos-stability-runge
{doc}`demos/stability-interpcheb`
::::

As a bonus, for Chebyshev nodes the barycentric weights are simple:

:::{math}
:label: weightcheb
  w_k = (-1)^k d_k, \qquad d_k =
  \begin{cases}
    1/2, & \text{if $k=0$ or $k=n$},\\
    1, & \text{otherwise}.
  \end{cases}
:::

## Spectral convergence

If we take $n\rightarrow \infty$ and use polynomial interpolation on Chebyshev nodes, the convergence rate is exponential in $n$. The following is typical of the results that can be proved.

::::{prf:theorem}
:label: theorem-spectral
Suppose $f(x)$ is analytic in an open real interval containing $[-1,1]$. Then there exist constants $C>0$ and $K>1$ such that
  
:::{math}
:label: spectral
\max_{x\in[-1,1]} | f(x) - p(x) | \le C K^{-n},
:::

where $p$ is the unique polynomial of degree $n$ or less defined by interpolation on $n+1$ Chebyshev 2nd kind points.
::::

The condition "$f$ is analytic" means that the Taylor series of $f$ converges to $f(x)$ in an open interval containing $[-1,1]$.[^analytic] A necessary condition of analyticity is that $f$ is infinitely differentiable.

[^analytic]: Alternatively, analyticity means that the function is extensible to one that is differentiable in the complex plane.

```{index} convergence rate; spectral
```
In some contexts we refer to {eq}`spectral` as linear convergence, but here it is typical to say that the rate is exponential, geometric, or {term}`spectral convergence`. One achieves  constant reduction factors in the error by constant increments of $n$. By contrast, algebraic convergence in the form $O(n^{-p})$ for some $p>0$ requires *multiplying* $n$ by a constant factor in order to reduce error by a constant factor. Graphically, spectral error is linear on a log–linear scale, while algebraic convergence is a straight line on a log–log scale.

## Exercises

1. ⌨ Revisit {prf:ref}`demos-stability-equispaced` and determine a good approximate value for the constant $K$ mentioned in the comments there. 

2. ⌨ For each case, compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=4,8,12,\ldots,60$. At each value of $n$, compute the infinity-norm error (that is, $\max |p(x)-f(x)|$ evaluated for at least 4000 values of $x$). Using a log–linear scale, plot the error as a function of $n$, then determine a good approximation to the constant $K$ in {eq}`spectral`. 

    **(a)** $f(x) = 1/(25x^2+1)\qquad$ 
    **(b)** $f(x) = \tanh(5 x+2)$

    **(c)** $f(x) = \cosh(\sin x)\qquad$
    **(d)** $f(x) = \sin(\cosh x)$

    ::::{only} solutions
    ```{code-block} matlab
    function ChebError

        function run_example
            xi = linspace(-1,1,801)';
            figure
            for n = [4 16 64]
                subplot(1,2,1)
                x = -cos(pi*(0:n)'/n);
                yi = chebinterp(xi,x,f(x));
                semilogy(xi,abs(yi-f(xi))), hold on
                subplot(1,2,2)
                x = linspace(-1,1,n+1)';
                p = polyinterp(x,f(x));
                yi = p(xi);
                semilogy(xi,abs(yi-f(xi))), hold on
            end
        end

    f = @(x) 1./(25*x.^2 + 1);
    run_example

    f = @tanh;
    run_example

    f = @(x) cosh(sinh(x));
    run_example

    f = @abs;
    run_example

    end


    function P = chebinterp(xi,x,y)

    n = length(y)-1;
    x = -cos((0:n)'*pi/n);
    wghts = (-1).^(0:n)';
    wghts([1 n+1]) = wghts([1 n+1])/2;
    P = nan(size(xi));
    for j = 1:numel(P)
        t = wghts ./ (x -xi(j)-eps);
        P(j) = y'*t / sum(t);
    end

    end
    ```
    ::::

    (problem-chebinterp)=
3. ⌨ Write a function `chebinterp(f,n)` that returns a function representing the polynomial interpolant of the input function `f` using $n+1$ Chebyshev second kind nodes over $[-1,1]$. You should use {eq}`weightcheb` to compute the barycentric weights directly, rather than using the method in {numref}`Function {number}<function-polyinterp>`. Test your function by revisiting {prf:ref}`demos-stability-runge` to use Chebyshev rather than equally spaced nodes. 

4. {prf:ref}`theorem-spectral` assumes that the function being approximated has infinitely many derivatives over $[-1,1]$. But now consider the family of functions $f_m(x)=|x|^m$. 

    **(a)** ✍ How many continuous derivatives over $[-1,1]$ does $f_m$ possess?

    **(b)** ⌨ Compute the polynomial interpolant using $n$ second-kind Chebyshev nodes in $[-1,1]$ for $n=10,20,30,\ldots,100$. At each value of $n$, compute the infinity-norm error (that is, $\max |p(x)-f_m(x)|$ evaluated for at least 41000 values of $x$). Using a single log–log (not log–linear!) graph, plot the error as a function of $n$ for all 6 values $m=1,3,5,7,9,11$.
    
    **(c)** ✍  Based on the results of parts (a) and (b), form a hypothesis about the asymptotic behavior of the error for fixed $m$ as $n\rightarrow \infty$. 

    ::::{only} solutions
    ```{code-block}
    %%
    %
    %% (a)
    xi = linspace(-1,1,1001)';
    clf
    for m = 1:2:11
        f = @(x) abs(x).^m;
        for n = 4:2:100
            x = -cos(pi*(0:n)'/n);
            p = polyinterp(x,f(x));
            err(m,n) = norm( f(xi) - p(xi), inf );
        end
        loglog(4:2:100,err(m,4:2:100),'.-'), hold on
    end
    axis tight

    %% (b)
    % Because the curves are straight lines, we should assume a relationship
    % of the form $E_m(n) = C n^p$ for different $C$ and $p$. Then $p$ is
    % the slope of the line on the log-log scales.  
    slope = [ ];
    for m = 1:2:11
        fit = polyfit( log(20:2:80), log(err(m,20:2:80)), 1 );
      slope  =  [ slope, fit(1) ];
    end
    slope
    ```
    ::::

    (problem-changeinterval)=
5. The Chebyshev points can be used when the interval of interpolation is $[a,b]$ rather than $[-1,1]$ by means of the change of variable

    :::{math}
    :label: changeinterval
      z = \psi(x) = a + (b-a)\frac{(x+1)}{2}.
    :::

    **(a)** ✍  Show that $\psi(-1) = a$, $\psi(1) = b$, and $\psi$ is strictly increasing on $[-1,1]$.

    **(b)** ✍ Invert the relation {eq}`changeinterval` to solve for $x$ in terms of $\psi^{-1}(z)$. 

    **(c)** ✍ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. Then a polynomial in $x$ can be used to interpolate the function values $f\bigl(\psi(t_i)\bigr)$. This in turn implies an interpolating function $\tilde{P}(z) =P\bigl(\psi^{-1}(z)\bigr)$. Show that $\tilde{P}$ is a polynomial in $z$. 
	
    **(d)** ⌨ Implement the idea of part (c) to plot a polynomial interpolant of $f(x) =\cosh(\sin x)$ over $[0,2\pi]$ using $n+1$ Chebyshev nodes with $n=40$. 

6. The Chebyshev points can be used for interpolation of functions defined on the entire real line by using the change of variable

    :::{math}
    :label: changeintervalinf
    z = \phi(x) = \frac{2x}{1-x^2},
    :::

    which maps the interval $(-1,1)$ in one-to-one fashion to the entire real line.

    **(a)** ✍ Find $\displaystyle \lim_{x\to 1^-} \phi(x)$ and $\displaystyle \lim_{x\to -1^+} \phi(x)$. 

    **(b)** ✍ Invert {eq}`changeintervalinf` to express $x=\phi^{-1}(z)$. (Be sure to enforce $-1\le x \le 1$.)

    **(c)** ⌨ Let $t_0,\ldots,t_n$ be standard second-kind Chebyshev points. These map to the $z$ variable as $\zeta_i=\phi(t_i)$ for all $i$. Suppose that $f(z)$ is a given function whose domain is the entire real line. Then the function values $y_i=f(\zeta_i)$ can be associated with the Chebyshev nodes $t_i$, leading to a polynomial interpolant $p(x)$. This in turn implies an interpolating function on the real line, defined as

    $$
    q(z)=p\bigl(\phi^{-1}(z)\bigr) = p(x).
    $$
    
    Implement this idea to plot an interpolant of $f(z)=\tanh(3z-2)$ using $n=30$. Your plot should show $q(z)$ evaluated at 1000 evenly spaced points in $[-5,5]$, with markers at the nodal values (those lying within the $[-5,5]$ window).

    ::::{only} solutions
    %% (b)
    % $$x = (\sqrt{1+z^2}-1)/z$$
    %% (c)
    t = -cos(pi*(0:30)'/30);
    zeta = 2*t./(1-t.^2);
    y = tanh(3*zeta-2);
    p = polyinterp(t,y);
    z = linspace(-5,5,1000)';
    x = (sqrt(z.^2+1)-1)./z;
    plot(z,p(x))
    hold on, plot(zeta,y,'o'), xlim([-5 5])
    ::::
