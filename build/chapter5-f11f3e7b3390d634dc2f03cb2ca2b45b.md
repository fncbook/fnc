---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 5

MATLAB implementations

## Functions


(function-hatfun-matlab)=
``````{dropdown} Hat function
:open:
```{literalinclude} FNC-matlab/hatfun.m
:language: matlab
:linenos: true
```
``````

(function-plinterp-matlab)=
``````{dropdown} Piecewise linear interpolation
:open:
```{literalinclude} FNC-matlab/plinterp.m
:language: matlab
:linenos: true
```
``````

(function-spinterp-matlab)=
``````{dropdown} Cubic spline interpolation
:open:
```{literalinclude} FNC-matlab/spinterp.m
:language: matlab
:linenos: true
```
``````

(function-fdweights-matlab)=
``````{dropdown} Fornberg's algorithm for finite difference weights
:open:
```{literalinclude} FNC-matlab/fdweights.m
:language: matlab
:linenos: true
```
``````

(function-trapezoid-matlab)=
``````{dropdown} Trapezoid formula for numerical integration
:open:
```{literalinclude} FNC-matlab/trapezoid.m
:language: matlab
:linenos: true
```
``````

(function-intadapt-matlab)=
``````{dropdown} Adaptive integration
:open:
```{literalinclude} FNC-matlab/intadapt.m
:language: matlab
:linenos: true
```
:::{admonition} About the code
:class: dropdown
The intended way for a user to call {numref}`Function {number} <function-intadapt>` is with only `f`, `a`, `b`, and `tol` provided. We then use default values on the other parameters to compute the function values at the endpoints, the interval's midpoint, and the function value at the midpoint. Recursive calls from within the function itself will provide all of that information, since it was already calculated along the way.
:::
``````

## Examples

```{code-cell}
:tags: [remove-cell]
cd /Users/driscoll/Dropbox/Mac/Documents/GitHub/fnc/matlab
FNC_init
```

### 5.1 @section-localapprox-interpolation

(demo-interpolation-global-matlab)=
``````{dropdown} @demo-interpolation-global
Here are some points that we could consider to be observations of an unknown function on $[-1,1]$.

```{code-cell}
n = 5;
t = linspace(-1,1,n+1)';  
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t,y)
```

```{index} ! Julia; fit
```

The polynomial interpolant, as computed using `polyfit`, looks very sensible. It's the kind of function you'd take home to meet your parents.

```{code-cell}
c = polyfit(t, y, n);     % polynomial coefficients
p = @(x) polyval(c, x);
hold on
fplot(p, [-1 1])
legend('data', 'interpolant', 'location', 'north');
```

But now consider a different set of points generated in almost exactly the same way.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1);
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t, y)
```

The points themselves are unremarkable. But take a look at what happens to the polynomial interpolant.

```{code-cell}
c = polyfit(t, y, n);     % polynomial coefficients
p = @(x) polyval(c, x);
hold on, fplot(p, [-1 1])
legend('data', 'interpolant', 'location', 'north');
```

Surely there must be functions that are more intuitively representative of those points!
``````

(demo-interpolation-pwise-matlab)=
``````{dropdown} @demo-interpolation-pwise
Let us recall the data from {numref}`Demo %s <demo-interpolation-global>`.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1);
y = t.^2 + t + 0.05 * sin(20 * t);
clf, scatter(t, y)
```

Here is an interpolant that is linear between each consecutive pair of nodes, using `interp1` from MATLAB.

```{code-cell}
x = linspace(-1, 1, 400)';
hold on, plot(x, interp1(t, y, x))
title('Piecewise linear interpolant') 
```

We may prefer a smoother interpolant that is piecewise cubic, generated using `Spline1D` from the `Dierckx` package.

```{code-cell}
cla
scatter(t, y)
plot(x, interp1(t, y, x, 'spline'))
title('Piecewise cubic interpolant')  
```
``````

(demo-interp-cond-matlab)=
``````{dropdown} @demo-interp-cond
In {numref}`Demo %s <demo-interpolation-global>` and {numref}`Demo %s <demo-interpolation-pwise>` we saw a big difference between polynomial interpolation and piecewise polynomial interpolation of some arbitrarily chosen data. The same effects can be seen clearly in the cardinal functions, which are closely tied to the condition numbers.

```{code-cell}
n = 18;
t = linspace(-1, 1, n+1)';
y = [zeros(9, 1); 1; zeros(n - 9, 1)];    % 10th cardinal function
clf, scatter(t, y)
hold on
x = linspace(-1, 1, 400)';
plot(x, interp1(t, y, x, 'spline'))
title('Piecewise cubic cardinal function') 
xlabel('x'), ylabel('p(x)') 
```

The piecewise cubic cardinal function is nowhere greater than one in absolute value. This happens to be true for all the cardinal functions, ensuring a good condition number for any interpolation with these functions. But the story for global polynomials is very different.

```{code-cell}
clf, scatter(t, y)
c = polyfit(t, y, n);
hold on, plot(x, polyval(c, x))
title('Polynomial cardinal function')
xlabel('x'), ylabel(('p(x)'));
```

From the figure we can see that the condition number for polynomial interpolation on these nodes is at least 500.
``````

### 5.2 @section-localapprox-pwlin

(demo-pwlin-hat-matlab)=
``````{dropdown} @demo-pwlin-hat
Let's define a set of four nodes (i.e., $n=3$ in our formulas).

```{index} ! Julia; annotate!
```

```{code-cell}
t = [0, 0.55, 0.7, 1];
```

We plot the hat functions $H_0,\ldots,H_3$.

```{code-cell}
clf
for k = 0:3
    subplot(4, 1, k+1)
    Hk = hatfun(t, k);
    fplot(Hk, [0, 1])
    hold on
    scatter(t, Hk(t))
    text(t(k+1), 0.6, sprintf("H_%d", k))
end
```
``````

(demo-pwlin-usage-matlab)=
``````{dropdown} @demo-pwlin-usage
We generate a piecewise linear interpolant of $f(x)=e^{\sin 7x}$.

```{code-cell}
f = @(x) exp(sin(7 * x));
clf
fplot(f, [0, 1], displayname="function")
xlabel("x");  ylabel(("y"));
```

First we sample the function to create the data.

```{code-cell}
t = [0, 0.075, 0.25, 0.55, 0.7, 1];    % nodes
y = f(t);                              % function values
```

Now we create a callable function that will evaluate the piecewise linear interpolant at any $x$, and then plot it.

```{code-cell}
p = plinterp(t, y);
hold on
fplot(p, [0, 1], displayname="interpolant")
scatter(t, y, displayname="values at nodes")
title("PL interpolation")
legend();
```
``````

(demo-pwlin-converge-matlab)=
``````{dropdown} @demo-pwlin-converge
We measure the convergence rate for piecewise linear interpolation of $e^{\sin 7x}$ over $x \in [0,1]$.

```{code-cell}
f = @(x) exp(sin(7 * x));
x = linspace(0, 1, 10001)';    % sample the difference at many points
n = round(10.^(1:0.25:3.5))';
maxerr = zeros(size(n));
for i = 1:length(n)
    t = (0:n(i)) / n(i);       % interpolation nodes
    p = plinterp(t, f(t));
    maxerr(i) = norm(f(x) - p(x), Inf);
end
disp(table(n(1:4:end), maxerr(1:4:end), variableNames=["n", "inf-norm error"]))
```

As predicted, a factor of 10 in $n$ produces a factor of 100 reduction in the error. In a convergence plot, it is traditional to have $h$ *decrease* from left to right, so we expect a straight line of slope $-2$ on a log-log plot.

```{code-cell}
clf
loglog(n, maxerr, "-o", displayname="error")
order2 = 0.5 * maxerr(end) * (n / n(end)) .^ (-2);
hold on
loglog(n, order2, "k--", displayname="O(n^{-2})")
xlabel("n");  ylabel("|| f-p ||_{\infty}")
title("Convergence of PL interpolation")
legend();
```
``````

### 5.3 @section-localapprox-splines 

(demo-splines-splines-matlab)=
``````{dropdown} @demo-splines-splines
For illustration, here is a spline interpolant using just a few nodes.

```{code-cell}
clf
f = @(x) exp(sin(7 * x));
fplot(f, [0, 1], displayname="function")

t = [0, 0.075, 0.25, 0.55, 0.7, 1];    % nodes
y = f(t);                              % values at nodes
hold on, scatter(t, y, displayname="values at nodes")

S = spinterp(t, y);
fplot(S, [0, 1], displayname="spline")

xlabel("x");  ylabel("y")
legend();
```

Now we look at the convergence rate as the number of nodes increases.

```{code-cell}
x = (0:10000)' / 1e4;              % sample the difference at many points
n = round(2 .^ (3:0.5:7))';        % numbers of nodes
maxerr = zeros(size(n));
for i = 1:length(n)
    t = (0:n(i))' / n(i);
    S = spinterp(t, f(t));
    err = f(x) - S(x);
    maxerr(i) = norm(err, Inf);
end
disp(table(n(1:2:end), maxerr(1:2:end), variableNames=["n", "inf-norm error"]))
```

Since we expect convergence that is $O(h^4)=O(n^{-4})$, we use a log-log graph of error and expect a straight line of slope $-4$.

```{code-cell}
clf
loglog(n, maxerr, "-o", displayname="error")
order4 = 0.5 * maxerr(end) * (n / n(end)) .^ (-4);
hold on
loglog(n, order4, "k--", displayname="O(n^{-4})")
xlabel("n");  ylabel("|| f-S ||_{\infty}")
title(("Convergence of spline interpolation"));
```
``````

### 5.4 @section-localapprox-finitediffs

(demo-finitediffs-fd1-matlab)=
``````{dropdown} @demo-finitediffs-fd1
If $f(x)=e^{\,\sin(x)}$, then $f'(0)=1$.

```{code-cell}
f = @(x) exp(sin(x));
```

Here are the first two centered differences from {numref}`table-FDcenter`.

```{code-cell}
h = 0.05;
format long
CD2 = (-f(-h) + f(h)) / (2*h)
CD4 = (f(-2*h) - 8*f(-h) + 8*f(h) - f(2*h)) / (12*h)
```

Here are the first two forward differences from {numref}`table-FDforward`.

```{code-cell}
FD1 = (-f(0) + f(h)) / h
FD2 = (-3*f(0) + 4*f(h) - f(2*h)) / (2*h)
```

Finally, here are the backward differences that come from reverse-negating the forward differences.

```{code-cell}
BD1 = (-f(-h) + f(0)) / h
BD2 = (f(-2*h) - 4*f(-h) + 3*f(0)) / (2*h)
```
``````

(demo-finitediffs-fd2-matlab)=
``````{dropdown} @demo-finitediffs-fd2
If $f(x)=e^{\,\sin(x)}$, then $f''(0)=1$.

```{code-cell}
f = @(x) exp(sin(x));
```

Here is a centered estimate given by {eq}`centerFD22`.

```{code-cell}
h = 0.05;
format long
CD2 = (f(-h) - 2*f(0) + f(h)) / h^2
```

For the same $h$, here are forward estimates given by {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
FD1 = (f(0) - 2*f(h) + f(2*h)) / h^2
FD2 = (2*f(0) - 5*f(h) + 4*f(2*h) - f(3*h)) / h^2
```

Finally, here are the backward estimates that come from reversing {eq}`forwardFD21` and {eq}`forwardFD22`.

```{code-cell}
BD1 = (f(-2*h) - 2*f(-h) + f(0)) / h^2
BD2 = (-f(-3*h) + 4*f(-2*h) - 5*f(-h) + 2*f(0)) / h^2
```
``````

(demo-finitediffs-fd-weights-matlab)=
``````{dropdown} @demo-finitediffs-fd-weights
We will estimate the derivative of $\cos(x^2)$ at $x=0.5$ using five nodes.

```{code-cell}
t = [0.35, 0.5, 0.57, 0.6, 0.75];    % nodes
f = @(x) cos(x.^2);
dfdx = @(x) -2 * x * sin(x^2);
exact_value = dfdx(0.5)
```

We have to shift the nodes so that the point of estimation for the derivative is at $x=0$. (To subtract a scalar from a vector, we must use the `.-` operator.)

```{code-cell}
format short
w = fdweights(t - 0.5, 1)
```

The finite-difference formula is a dot product (i.e., inner product) between the vector of weights and the vector of function values at the nodes.

```{code-cell}
fd_value = w * f(t)'
```

We can reproduce the weights in the finite-difference tables by using equally spaced nodes with $h=1$. For example, here is a one-sided formula at four nodes.

```{code-cell}
fdweights(0:3, 1)
```

``````

### 5.5 @section-localapprox-fd-converge

(demo-fdconverge-order12-matlab)=
``````{dropdown} @demo-fdconverge-order12
Let's observe the convergence of the formulas in {numref}`Example {number} <example-fd-converge-FD11>` and {numref}`Example {number} <example-fd-converge-FD12>`, applied to the function $\sin(e^{x+1})$ at $x=0$.

```{code-cell}
f = @(x) sin(exp(x + 1));
exact_value = exp(1) * cos(exp(1))
```

We'll compute the formulas in parallel for a sequence of $h$ values.

```{code-cell}
h = 5 ./ 10.^(1:6)';
FD1 = zeros(size(h));
FD2 = zeros(size(h));
for i = 1:length(h)
    h_i = h(i);
    FD1(i) = (f(h_i) - f(0)    ) / h_i;
    FD2(i) = (f(h_i) - f(-h_i)) / (2*h_i);
end
disp(table(h, FD1, FD2))
```

All that's easy to see from this table is that FD2 appears to converge to the same result as FD1, but more rapidly. A table of errors is more informative.

```{code-cell}
err1 = abs(exact_value - FD1);
err2 = abs(exact_value - FD2);
disp(table(h, err1, err2, variableNames=["h", "error in FD1", "error in FD2"]))
```

In each row, $h$ is decreased by a factor of 10, so that the error is reduced by a factor of 10 in the first-order method and 100 in the second-order method.

A graphical comparison can be useful as well. On a log-log scale, the error should (as $h\to 0$) be a straight line whose slope is the order of accuracy. However, it's conventional in convergence plots to show $h$ _decreasing_ from left to right, which negates the slopes.

```{code-cell}
clf
loglog(h, abs([err1 err2]), "o-")
set(gca, "xdir", "reverse")
order1 = 0.1 * err1(end) * (h / h(end)) .^ 1;
order2 = 0.1 * err2(end) * (h / h(end)) .^ 2;
hold on
loglog(h, order1, "--", h, order2, "--")
xlabel("h");  ylabel("error")
title("Convergence of finite differences")
legend("FD1", "FD2", "O(h)", "O(h^2)");
```
``````

(demo-fdconverge-round-matlab)=
``````{dropdown} @demo-fdconverge-round
Let $f(x)=e^{-1.3x}$. We apply finite-difference formulas of first, second, and fourth order to estimate $f'(0)=-1.3$.

```{code-cell}
f = @(x) exp(-1.3 * x);
exact = -1.3;

h = 10 .^ (-(1:12))';
FD = zeros(length(h), 3);
for i = 1:length(h)
    h_i = h(i);
    nodes = h_i * (-2:2);
    vals = f(nodes);
    FD(i, 1) = dot([0      0 -1   1    0] / h_i, vals);
    FD(i, 2) = dot([0    -1/2 0 1/2    0] / h_i, vals);
    FD(i, 3) = dot([1/12 -2/3 0 2/3 -1/12] / h_i, vals);
end
format long
disp(table(h, FD(:, 1), FD(:, 2), FD(:, 3), variableNames=["h", "FD1", "FD2", "FD4"]))
```

They all seem to be converging to $-1.3$. The convergence plot reveals some interesting structure to the errors, though.

```{code-cell}
err = abs(FD - exact);
clf
loglog(h, err, "o-")
set(gca, "xdir", "reverse")
order1 = 0.1 * err(end, 1) * (h / h(end)) .^ (-1);
hold on
loglog(h, order1, "k--")
xlabel("h");  ylabel("error")
title("FD error with roundoff")
legend("FD1", "FD2", "FD4", "O(1/h)", "location", "northeast");
```

Again the graph is made so that $h$ decreases from left to right. The errors are dominated at first by truncation error, which decreases most rapidly for the fourth-order formula. However, increasing roundoff error eventually equals and then dominates the truncation error as $h$ continues to decrease. As the order of accuracy increases, the crossover point moves to the left (greater efficiency) and down (greater accuracy).
``````

### 5.6 @section-localapprox-integration
(demo-int-antideriv-matlab)=
``````{dropdown} @demo-int-antideriv
The antiderivative of $e^x$ is, of course, itself. That makes evaluation of $\int_0^1 e^x\,dx$ by the Fundamental Theorem trivial.

```{code-cell}
format long
exact = exp(1) - 1
```

```{index} ! MATLAB; integral
```

MATLAB has numerical integrator `integral` that estimates the value without finding the antiderivative first. As you can see here, it can be as accurate as floating-point precision allows.

```{code-cell}
integral(@(x) exp(x), 0, 1)
```

The numerical approach is also far more robust. For example, $e^{\,\sin x}$ has no useful antiderivative. But numerically, it's no more difficult.

```{code-cell}
integral(@(x) exp(sin(x)), 0, 1)
```

When you look at the graphs of these functions, what's remarkable is that one of these areas is basic calculus while the other is almost impenetrable analytically. From a numerical standpoint, they are practically the same problem.

```{code-cell}
:tags: hide-input
x = linspace(0, 1, 201)';
subplot(2,1,1), fill([x; 1; 0], [exp(x); 0;0 ], [1, 0.9, 0.9])
title('exp(x)')
ylabel('f(x)')
subplot(2, 1, 2), fill([x; 1; 0], [exp(sin(x)); 0; 0], [1, 0.9, 0.9])
title('exp(sin(x))')
xlabel('x'), ylabel(('f(x)'));
```
``````

(demo-int-trap-matlab)=
``````{dropdown} @demo-int-trap
We will approximate the integral of the function $f(x)=e^{\sin 7x}$ over the interval $[0,2]$.

```{code-cell}
f = @(x) exp(sin(7 * x));
a = 0;  b = 2;
```

In lieu of the exact value, we use the `integral` function to find an accurate result.

```{code-cell}
I = integral(f, a, b, abstol=1e-14, reltol=1e-14);
fprintf("Integral = %.15f", I)
```

Here is the trapezoid result at $n=40$, and its error.

```{code-cell}
T = trapezoid(f, a, b, 40);
fprintf("Trapezoid error = %.2e", I - T)
```

In order to check the order of accuracy, we increase $n$ by orders of magnitude and observe how the error decreases.

```{code-cell}
n = 10 .^ (1:5)';
err = zeros(size(n));
for i = 1:length(n)
    T = trapezoid(f, a, b, n(i));
    err(i) = I - T;
end
disp(table(n, err, variableNames=["n", "Trapezoid error"]))
```

Each increase by a factor of 10 in $n$ cuts the error by a factor of about 100, which is consistent with second-order convergence. Another check is that a log-log graph should give a line of slope $-2$ as $n\to\infty$.

```{code-cell}
clf
loglog(n, abs(err), "-o", displayname="trapezoid")
hold on
loglog(n, 0.1 * abs(err(end)) * (n / n(end)).^(-2), "k--", displayname="O(n^{-2})")
xlabel("n");  ylabel("error")
title("Convergence of trapezoidal integration")
legend();
```
``````

(demo-int-extrap-matlab)=
``````{dropdown} @demo-int-extrap
We estimate $\displaystyle\int_0^2 x^2 e^{-2x}\, dx$ using extrapolation. First we use `quadgk` to get an accurate value.

```{code-cell}
f = @(x) x.^2 .* exp(-2 * x);
a = 0;  b = 2;
format long
I = integral(f, a, b, abstol=1e-14, reltol=1e-14)
```

We start with the trapezoid formula on $n=N$ nodes.

```{code-cell}
N = 20;       % the coarsest formula
n = N;  h = (b - a) / n;
t = h * (0:n)';
y = f(t);
```

We can now apply weights to get the estimate $T_f(N)$.

```{code-cell}
T = h * ( sum(y(2:n)) + y(1) / 2 + y(n+1) / 2 )
```

Now we double to $n=2N$, but we only need to evaluate $f$ at every other interior node and apply {eq}`nc-doubling`.

```{code-cell}
n = 2*n;  h = h / 2;
t = h * (0:n)';
T(2) = T(1) / 2 + h * sum( f(t(2:2:n)) )
```

We can repeat the same code to double $n$ again.

```{code-cell}
n = 2*n;  h = h / 2;
t = h * (0:n)';
T(3) = T(2) / 2 + h * sum( f(t(2:2:n)) )
```

Let us now do the first level of extrapolation to get results from Simpson's formula. We combine the elements `T[i]` and `T[i+1]` the same way for $i=1$ and $i=2$.

```{code-cell}
S = (4 * T(2:3) - T(1:2)) / 3
```

With the two Simpson values $S_f(N)$ and $S_f(2N)$ in hand, we can do one more level of extrapolation to get a sixth-order accurate result.

```{code-cell}
R = (16*S(2) - S(1)) / 15
```

We can make a triangular table of the errors:

```{code-cell}
err2 = T(:) - I;
err4 = [NaN; S(:) - I];
err6 = [NaN; NaN; R - I];
format short e
disp(table(err2, err4, err6, variablenames=["order 2", "order 4", "order 6"]))
```

If we consider the computational time to be dominated by evaluations of $f$, then we have obtained a result with about twice as many accurate digits as the best trapezoid result, at virtually no extra cost.
``````

### 5.7 @section-localapprox-adaptive

(demo-adapt-motive-matlab)=
``````{dropdown} @demo-adapt-motive
This function gets increasingly oscillatory as $x$ increases.

```{code-cell}
f = @(x) (x + 1).^2 .* cos((2 * x + 1) ./ (x - 4.3));
clf
fplot(f, [0, 4], 2000)
xlabel('x'), ylabel(('f(x)'));
```

Accordingly, the trapezoid rule is more accurate on the left half of this interval than on the right half.

```{code-cell}
left_val = integral(f, 0, 2, abstol=1e-14, reltol=1e-14);
right_val = integral(f, 2, 4, abstol=1e-14, reltol=1e-14);

n = round(50 * 2 .^ (0:3)');
err = zeros(length(n), 2);
for i = 1:length(n)
    T = trapezoid(f, 0, 2, n(i));
    err(i, 1) = T - left_val;
    T = trapezoid(f, 2, 4, n(i));
    err(i, 2) = T - right_val;
end
disp(table(n, err(:, 1), err(:, 2), variableNames=["n", "left error", "right error"]))
```

Both the picture and the numerical results suggest that more nodes should be used on the right half of the interval than on the left half.
``````

(demo-adapt-usage-matlab)=
``````{dropdown} @demo-adapt-usage
We'll integrate the function from {numref}`Demo %s <demo-adapt-motive>`.

```{code-cell}
f = @(x) (x + 1).^2 .* cos((2 * x + 1) ./ (x - 4.3));
```

We perform the integration and show the nodes selected underneath the curve.

```{code-cell}
[Q, t] = intadapt(f, 0, 4, 0.001);
clf, fplot(f, [0, 4], 2000)
hold on
stem(t, f(t), '.-')
title('Adaptive node selection')
xlabel('x'), ylabel('f(x)')
fprintf("number of nodes = %d", length(t))
```

The error turns out to be a bit more than we requested. It's only an estimate, not a guarantee.

```{code-cell}
I = integral(f, 0, 4, abstol=1e-14, reltol=1e-14);    % 'exact' value
fprintf("error = %.2e", abs(Q - I))
```

Let's see how the number of integrand evaluations and the error vary with the requested tolerance.

```{code-cell}
tol = 1 ./ 10.^(4:14)';
err = zeros(size(tol));
n = zeros(size(tol));
for i = 1:length(tol)
    [A, t] = intadapt(f, 0, 4, tol(i));
    err(i) =  I - A;
    n(i) = length(t);
end
disp(table(tol, err, n, variableNames=["tolerance", "error", "number of nodes"]))
```

As you can see, even though the errors are not smaller than the estimates, the two columns decrease in tandem. If we consider now the convergence not in $h$, which is poorly defined now, but in the number of nodes actually chosen, we come close to the fourth-order accuracy of the underlying Simpson scheme.

```{code-cell}
clf
loglog(n, abs(err), "-o", displayname="results")
xlabel("number of nodes"), ylabel("error")
title("Convergence of adaptive integration")
order4 = 0.1 * abs(err(end)) * (n / n(end)).^(-4);
hold on
loglog(n, order4, "k--", displayname="O(n^{-4})")
legend();
```
``````