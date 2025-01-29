---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
    headings: false
---
# Chapter 4

## Functions

(function-newton-matlab)=
``````{dropdown} Newton's method
:open:
```{literalinclude} FNC-matlab/newton.m
:language: matlab
:linenos: true
```
``````

(function-secant-matlab)=
``````{dropdown} Secant method
:open:
```{literalinclude} FNC-matlab/secant.m
:language: matlab
:linenos: true
```
``````

(function-newtonsys-matlab)=
``````{dropdown} Newton's method for systems
:open:
```{literalinclude} FNC-matlab/newtonsys.m
:language: matlab
:linenos: true
```
``````

(function-fdjac-matlab)=
``````{dropdown} Finite differences for Jacobian
:open:
```{literalinclude} FNC-matlab/fdjac.m
:language: matlab
:linenos: true
```
``````

(function-levenberg-matlab)=
``````{dropdown} Levenberg's method
:open:
```{literalinclude} FNC-matlab/levenberg.m
:language: matlab
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-cell]
cd  /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init;
```

### 4.1 @section-nonlineqn-rootproblem

(demo-rootproblem-bessel-matlab)=
``````{dropdown} @demo-rootproblem-bessel
:open:

```{code-cell}
J3 = @(x) besselj(3,x);
fplot(J3, [0, 20])
grid on
xlabel('x'), ylabel('J_3(x)')  
title('Bessel function') 
```
From the graph we see roots near 6, 10, 13, 16, and 19. We use `nlsolve` from the `NLsolve` package to find these roots accurately. It uses vector variables, so we have to code accordingly.
```{tip}
:class: dropdown
Type `\omega` followed by <kbd>Tab</kbd> to get the character `ω`.
The argument `ftol=1e-14` below is called a **keyword argument**. Here it sets a goal for the maximum value of $|f(x)|$.
```

```{code-cell}
omega = [];
for guess = [6, 10, 13, 16, 19]
    omega = [omega; fzero(J3, guess)];
end
omega
```

```{code-cell}
table(omega, J3(omega), 'VariableNames', {'root estimate', 'function value'})
```

```{code-cell}
hold on
scatter(omega, J3(omega))
title('Bessel roots')    
```

If instead we seek values at which $J_3(x)=0.2$, then we must find roots of the function $J_3(x)-0.2$.

```{code-cell}
omega = [];
for guess = [3, 6, 10, 13]
    f = @(x) J3(x) - 0.2;
    omega = [omega; fzero(f, guess)];
end
scatter(omega, J3(omega), '<')
```
``````

(demo-roots-cond-matlab)=
``````{dropdown} @demo-roots-cond
:open:
Consider first the function

```{code-cell}
f  = @(x) (x - 1) .* (x - 2);
```

At the root $r=1$, we have $f'(r)=-1$. If the values of $f$ were perturbed at every point by a small amount of noise, we can imagine finding the root of the function drawn with a thick ribbon, giving a range of potential roots.


```{code-cell}
clf
interval = [0.8, 1.2];
fplot(f, interval)
grid on, hold on
fplot(@(x) f(x) + 0.02, interval, 'k')
fplot(@(x) f(x) - 0.02, interval, 'k')
axis equal,  yticks([0])
ylim([-0.1, 0.1])
xlabel('x'), ylabel('f(x)')
title('Well-conditioned root')
```

The possible values for a perturbed root all lie within the interval where the ribbon intersects the $x$-axis. The width of that zone is about the same as the vertical thickness of the ribbon.

By contrast, consider the function

```{code-cell}
f = @(x) (x - 1) .* (x - 1.01);
```

Now $f'(1)=-0.01$, and the graph of $f$ will be much shallower near $x=1$. Look at the effect this has on our thick rendering:

```{code-cell}
axis(axis), cla
fplot(f, interval)
fplot(@(x) f(x) + 0.02, interval, 'k')
fplot(@(x) f(x) - 0.02, interval, 'k')
ylim([-0.1, 0.1]), yticks([0])
title('Poorly conditioned root') 
```

The vertical displacements in this picture are exactly the same as before. But the potential _horizontal_ displacement of the root is much wider. In fact, if we perturb the function entirely upward by the amount drawn here, the root disappears!
``````

### 4.2 @section-nonlineqn-fixed-point
(demo-fp-spiral-matlab)=
``````{dropdown} @demo-fp-spiral
:open:
Let's convert the roots of a quadratic polynomial $f(x)$ to a fixed point problem.

```{code-cell}
f = @(x) x.^2 - 4*x + 3.5;
r = roots([1, -4, 3.5])
```

We define $g(x)=x-p(x)$. 

```{code-cell}
g = @(x) x - f(x);
```

Intersections of $y=g(x)$ with the line $y=x$ are fixed points of $g$ and thus roots of $f$. (Only one is shown in the chosen plot range.)

```{code-cell}
clf
fplot(g, [2, 3])
hold on,  plot([2, 3], [2, 3], 'k')
title('Finding a fixed point'),  axis equal  
xlabel('x'),  ylabel('y')  
```

If we evaluate $g(2.1)$, we get a value of almost 2.6, so this is not a fixed point.

```{code-cell}
x = 2.1;
y = g(x)
```

However, $y=g(x)$ is considerably closer to the fixed point at around 2.7 than $x$ is. Suppose then that we adopt $y$ as our new $x$ value. Changing the $x$ coordinate in this way is the same as following a horizontal line over to the graph of $y=x$.

```{code-cell}
plot([x, y], [y, y], '-')
x = y;
```

Now we can compute a new value for $y$. We leave $x$ alone here, so we travel along a vertical line to the graph of $g$.

```{code-cell}
y = g(x)
plot([x, x],[x, y], '-')
```

You see that we are in a position to repeat these steps as often as we like. Let's apply them a few times and see the result.

```{code-cell}
for k = 1:5
    plot([x, y], [y, y], '-')
    x = y;       % y --> new x
    y = g(x);    % g(x) --> new y
    plot([x, x], [x, y], '-')  
end
```

The process spirals in beautifully toward the fixed point we seek. Our last estimate has almost 4 accurate digits.

```{code-cell} 
abs(y - r(1)) / r(1)
```

Now let's try to find the other fixed point $\approx 1.29$ in the same way. We'll use 1.3 as a starting approximation.

```{code-cell}
cla
fplot(g, [1, 2])
hold on, plot([1, 2], [1, 2], 'k')
ylim([1, 2])
x = 1.3;  y = g(x);
for k = 1:5
    plot([x, y], [y, y], '-'),  
    x = y;       % y --> new x
    y = g(x);    % g(x) --> new y
    plot([x, x], [x, y], '-') 
end
title('No convergence') 
```

This time, the iteration is pushing us _away from_ the correct answer.
``````

(demo-fp-converge-matlab)=
``````{dropdown} @demo-fp-converge
:open:
We revisit {numref}`Demo %s <demo-fp-spiral>` and investigate the observed convergence more closely. Recall that above we calculated $g'(p)\approx-0.42$ at the convergent fixed point.

```{code-cell}
f = @(x) x.^2 - 4*x + 3.5;
r = roots([1, -4, 3.5]);
```

Here is the fixed-point iteration. This time we keep track of the whole sequence of approximations.

:::{index} Julia; push!
:::

```{code-cell}
g = @(x) x - f(x);
x = 2.1; 
for k = 1:12
    x(k+1) = g(x(k));
end
```

It's illuminating to construct and plot the sequence of errors.

```{code-cell}
err = abs(x - r(1));
clf
semilogy(err, 'o-'), axis tight
xlabel('iteration'),  ylabel('error')
title('Convergence of fixed-point iteration')
```

It's quite clear that the convergence quickly settles into a linear rate. We could estimate this rate by doing a least-squares fit to a straight line. Keep in mind that the values for small $k$ should be left out of the computation, as they don't represent the linear trend.

```{code-cell}
y = log(err(5:12));
p = polyfit(5:12, y, 1);
```

We can exponentiate the slope to get the convergence constant $\sigma$.

```{code-cell}
sigma = exp(p(1))
```

The error should therefore decrease by a factor of $\sigma$ at each iteration. We can check this easily from the observed data.

```{code-cell}
err(9:12) ./ err(8:11)
```

The methods for finding $\sigma$ agree well.
``````
### 4.3 @section-nonlineqn-newton
(demo-newton-line-matlab)=
``````{dropdown} @demo-newton-line
:open:

Suppose we want to find a root of the function

```{code-cell}
f = @(x) x .* exp(x) - 2;
clf, fplot(f, [0, 1.5])
xlabel('x'), ylabel('y')    
set(gca, 'ygrid', 'on')  
title('Objective function')    

```

From the graph, it is clear that there is a root near $x=1$. So we call that our initial guess, $x_1$.

```{code-cell}
x1 = 1;
y1 = f(x1)
hold on, scatter(x1, y1, 'k')
```

Next, we can compute the tangent line at the point $\bigl(x_1,f(x_1)\bigr)$, using the derivative.

```{code-cell}
df_dx = @(x) exp(x) .* (x + 1);
slope1 = df_dx(x1);
tangent1 = @(x) y1 + slope1 * (x - x1);
axis(axis)
fplot(tangent1, [0, 1.5], 'k--')
title('Function and tangent line')    

```

In lieu of finding the root of $f$ itself, we settle for finding the root of the tangent line approximation, which is trivial. Call this $x_2$, our next approximation to the root.

```{code-cell}
x2 = x1 - y1 / slope1
scatter(x2, 0, 'r')
title('Root of the tangent')    
```

```{code-cell}
y2 = f(x2)
```

The residual (i.e., value of $f$) is smaller than before, but not zero. So we repeat the process with a new tangent line based on the latest point on the curve.

```{code-cell}
cla,  axis auto
fplot(f, [0.83, 0.88])
scatter(x2, y2, 'k')
slope2 = df_dx(x2);
tangent2 = @(x) y2 + slope2 * (x - x2);
axis(axis)
fplot(tangent2, [0.8, 0.9], 'k--')
x3 = x2 - y2 / slope2;
scatter(x3, 0, 'r')
title('Next iteration')    
```

```{code-cell}
y3 = f(x3)
```

Judging by the residual, we appear to be getting closer to the true root each time.
``````

(demo-newton-converge-matlab)=
``````{dropdown} @demo-newton-converge
:open:
We again look at finding a solution of $x e^x=2$ near $x=1$. To apply Newton's method, we need to calculate values of both the residual function $f$ and its derivative.

```{code-cell}
f = @(x) x.*exp(x) - 2;
df_dx = @(x) exp(x).*(x+1);
```

We don't know the exact root, so we use `nlsolve` to determine a proxy for it.

```{code-cell}
format long,  r = fzero(f,1)
```

We use $x_1=1$ as a starting guess and apply the iteration in a loop, storing the sequence of iterates in a vector.

```{code-cell}
x = 1;
for k = 1:6
    x(k+1) = x(k) - f(x(k)) / df_dx(x(k));
end
x
```

Here is the sequence of errors.

```{code-cell}
format short e
err = x' - r
```

The exponents in the scientific notation definitely suggest a squaring sequence. We can check the evolution of the ratio in {eq}`quadratictest`.

```{code-cell}
format short
logerr = log(abs(err))
```

The clear convergence to 2 above constitutes good evidence of quadratic convergence.
``````

(demo-newton-usage-matlab)=
``````{dropdown} @demo-newton-usage
:open:
```{index} ! Julia; enumerate
```

Suppose we want to evaluate the inverse of the function $h(x)=e^x-x$. This means solving $y=e^x-x$ for $x$ when $y$ is given, which has no elementary form. If a value of $y$ is given numerically, though, we simply have a rootfinding problem for $f(x)=e^x-x-y$.
```{tip}
:class: dropdown
When a function is created, it can refer to any variables in scope at that moment. Those values are locked in to the definition, which is called a _closure_. If the enclosed variables change values later, the function still uses the values it was created with.
```

```{code-cell}
h = @(x) exp(x) - x;
dh_dx = @(x) exp(x) - 1;
y_ = linspace(h(0), h(2), 200);
x_ = zeros(size(y_));
for i = 1:length(y_)
    f = @(x) h(x) - y_(i);
    df_dx = @(x) dh_dx(x);
    x = newton(f, df_dx, 1);  x_(i) = x(end);
end
```

```{code-cell}
:tags: [hide-input]
clf, fplot(h, [0, 2])
hold on, axis equal
plot(y_, x_)
plot([0, max(y_)], [0, max(y_)], 'k--')
xlabel('x'), ylabel('y')
legend('h(x)', 'inverse', 'y=x');
```
``````

### 4.4 @section-nonlineqn-secant

(demo-secant-line-matlab)=
``````{dropdown} @demo-secant-line
:open:


We return to finding a root of the equation $x e^x=2$.

```{code-cell}
f = @(x) x .* exp(x) - 2;
clf, fplot(f, [0.25, 1.25])
set(gca, 'ygrid', 'on')  
xlabel('x'), ylabel('y')    
title('Objective function')    
```

From the graph, it's clear that there is a root near $x=1$. To be more precise, there is a root in the interval $[0.5,1]$. So let us take the endpoints of that interval as _two_ initial approximations.

```{code-cell}
x1 = 1;    y1 = f(x1);
x2 = 0.5;  y2 = f(x2);
hold on, scatter([x1, x2], [y1, y2])
title('Two initial values')    
```

Instead of constructing the tangent line by evaluating the derivative, we can construct a linear model function by drawing the line between the two points $\bigl(x_1,f(x_1)\bigr)$ and $\bigl(x_2,f(x_2)\bigr)$. This is called a _secant line_.

```{code-cell}
slope2 = (y2 - y1) / (x2 - x1);
secant2 = @(x) y2 + slope2 * (x - x2);
```

As before, the next root estimate in the iteration is the root of this linear model.

```{code-cell}
fplot(secant2,[0.25, 1.25],'k--')
x3 = x2 - y2 / slope2;
y3 = f(x3)
scatter(x3, 0)
title('Next value') 
```

For the next linear model, we use the line through the two most recent points. The next iterate is the root of that secant line, and so on.

```{code-cell}
slope2 = (y3 - y2) / (x3 - x2);
x4 = x3 - y3 / slope2;
y4 = f(x4)
```
``````

(demo-secant-converge-matlab)=
``````{dropdown} @demo-secant-converge
:open:
We check the convergence of the secant method from {numref}`Demo %s <demo-secant-line>`. 

```{code-cell}
f = @(x) x .* exp(x) - 2;
x = secant(f, 1, 0.5);
```

We don't know the exact root, so we use `fzero` to get a good proxy.

```{code-cell}
r = fzero(f, 1);
```

Here is the sequence of errors.

```{code-cell}
format short e
err = r - x(1:end-1)'
```

It's not easy to see the convergence rate by staring at these numbers. We can use {eq}`superlinear-rate` to try to expose the superlinear convergence rate.

```{code-cell}
logerr = log(abs(err));
ratios = logerr(2:end) ./ logerr(1:end-1)
```

As expected, this settles in at around 1.618.
``````

(demo-secant-iqi-matlab)=
``````{dropdown} @demo-secant-iqi
:open:
Here we look for a root of $x+\cos(10x)$ that is close to 1.

```{code-cell}
f = @(x) x + cos(10 * x);
interval = [0.5, 1.5];
clf, fplot(f, interval)
set(gca, 'ygrid', 'on'), axis(axis)   
title('Objective function')    
xlabel('x'), ylabel('y')    
r = fzero(f, 1)
```

We choose three values to get the iteration started.

```{code-cell}
x = [0.8, 1.2, 1]';
y = f(x);
hold on, scatter(x, y)
title('Three initial points')    
```

If we were using forward interpolation, we would ask for the polynomial interpolant of $y$ as a function of $x$. But that parabola has no real roots.

```{code-cell}
c = polyfit(x, y, 2);    % coefficients of interpolant
q = @(x) polyval(c, x);
fplot(q, interval, '--')
title('Parabola model')     
```

To do inverse interpolation, we swap the roles of $x$ and $y$ in the interpolation.

```{tip}
:class: dropdown
By giving two functions in the `fplot` call, we get the parametric plot $(q(y),y)$ as a function of $y$.
```

```{code-cell}
cla, fplot(f, interval)
scatter(x, y)     
c = polyfit(y, x, 2);    % coefficients of interpolating polynomial
q = @(y) polyval(c, y);
fplot(q, @(y) y, ylim,'--')    % plot x=q(y), y=y
title('Sideways parabola')    
```

We seek the value of $x$ that makes $y$ zero. This means evaluating $q$ at zero.

```{code-cell}
x = [x; q(0)];
y = [y; f(x(end))]
```

We repeat the process a few more times.

```{code-cell}
for k = 4:8
    c = polyfit(y(k-2:k), x(k-2:k), 2);
    x(k+1) = polyval(c, 0);
    y(k+1) = f(x(k+1));
end
disp('final residual:')
y(end)
```

Here is the sequence of errors.

```{code-cell}
format short e
err = x - r
```

The convergence is probably superlinear:

```{code-cell}
logerr = log(abs(err));
ratios = logerr(2:end) ./ logerr(1:end-1)
```

``````
### 4.5 @section-nonlineqn-newtonsys
(demo-newtonsys-converge-matlab)=
``````{dropdown} @demo-newtonsys-converge
:open:
A system of nonlinear equations is defined by its residual and Jacobian.
```{tip}
:class: dropdown
This function needs to be defined within a script file or in a file of its own with the `.m` extension.
```

```{literalinclude} f45_nlsystem.m
:language: matlab
```

Since our system function is defined in an external file here, we need to use `@` in order to reference it as a function argument. 

```{code-cell}
nlsystem = @f45_nlsystem;
x1 = [0; 0; 0];    % column vector!
x = newtonsys(nlsystem, x1);
num_iter = size(x, 2)
```

Let's compute the residual of the last result in order to check the quality.

```{code-cell}
r = x(:, end)
back_err = norm(nlsystem(r))
```

We take the sequence errors in the first component of the solution, applying the log so that we can look at the exponents.

```{code-cell}
log10( abs(x(1, 1:end-1) - r(1)) )'
```

This sequence looks to be nearly doubling at each iteration, which is a good sign of quadratic convergence.
``````

### 4.6 @section-nonlineqn-quasinewton
(demo-quasi-levenberg-matlab)=
``````{dropdown} @demo-quasi-levenberg
:open:

To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.
```{tip}
:class: dropdown
A rule of thumb is that if you use a function as an input argument for another function, there needs to be an `@` involved once: either for an anonymous definition or to reference a function defined elsewhere. 
```{literalinclude} f45_nlsystem.m
:language: matlab
```
In all other respects usage is the same as for the `newtonsys` function.
```{code-cell}
f = @f46_nlsystem;
x1 = [0; 0; 0];   
x = levenberg(f, x1);
```
It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).
```{code-cell}
r = x(:, end)
backward_err = norm(f(r))
```
Looking at the convergence of the first component, we find a rate between linear and quadratic, like with the secant method.
```{code-cell}
log10( abs(x(1, 1:end-1) - r(1)) )'
```
``````
### 4.7 @section-nonlineqn-nlsq
(demo-nlsq-converge-matlab)=
``````{dropdown} @demo-nlsq-converge
:open:
We will observe the convergence of {numref}`Function {number} <function-levenberg>` for different levels of the minimum least-squares residual. We start with a function mapping from $\real^2$ into $\real^3$, and a point that will be near the optimum.
```{code-cell}
g = @(x) [sin(x(1) + x(2)); cos(x(1) - x(2)); exp(x(1) - x(2))];
p = [1; 1];
```
```{index} ! Julia; @sprintf
```

The function $\mathbf{g}(\mathbf{x}) - \mathbf{g}(\mathbf{p})$ obviously has a zero residual at $\mathbf{p}$. We'll make different perturbations of that function in order to create nonzero residuals.

```{tip}
:class: dropdown
`@sprintf` is a way to format numerical values as strings, patterned after the C function `printf`.
```

```{code-cell}
clf
labels = [];
for R = [1e-3, 1e-2, 1e-1]
    % Define the perturbed function.
    f = @(x) g(x) - g(p) + R * [-1; 1; -1] / sqrt(3)
    x = levenberg(f, [0; 0]);
    r = x(:, end);
    err = abs(x(1, 1:end-1) - r(1));
    normres = norm(f(r));
    semilogy(err), hold on
    labels = [labels; sprintf("R=%.2g", normres)];
end
xlabel("iteration"), ylabel("error")
legend(labels);
```

In the least perturbed case, where the minimized residual is less than $10^{-3}$, the convergence is plausibly quadratic. At the next level up, the convergence starts similarly but suddenly stagnates for a long time. In the most perturbed case, the quadratic phase is nearly gone and the overall shape looks linear.
``````

(demo-nlsq-MM-matlab)=
``````{dropdown} @demo-nlsq-MM
:open:
```{code-cell}
m = 25; V = 2; Km = 0.5;
s = linspace(0.05, 6, m)';
w = V * s ./ (Km + s);                      % exactly on the curve
w = w + 0.15 * cos(2 * exp(s / 16) .* s);   % noise added
clf, fplot(@(s) V * s ./ (Km + s), [0, 6], '--')
hold on, scatter(s, w)
xlabel('concentration'), ylabel('reaction rate')    
labels = ["ideal", "noisy data"];    
legend(labels, 'location', 'east');
```

The idea is to pretend that we know nothing of the origins of this data and use nonlinear least squares to find the parameters in the theoretical model function $v(s)$. In {eq}`nlsq-misfit`, the $s$ variable plays the role of $t$, and $v$ plays the role of $g$. In the Jacobian, the derivatives are with respect to the parameters in $\mathbf{x}$.

```{literalinclude} f47_misfit.m
:language: matlab
```

The misfit function above has to know the parameters `x` that are being optimized as well as the data `s` and `w` that remain fixed. We use a closure to pass the data values along.

```{code-cell}
f = @(x) f47_misfit(x, s, w);
```

Now we have a function that accepts a single 2-vector input and returns a 25-vector output. We can pass this function to `levenberg` to find the best-fit parameters.

```{code-cell}
x1 = [1; 0.75];
x = newtonsys(f, x1);
V = x(1, end),  Km = x(2, end)     % final values
model = @(s) V * s ./ (Km + s);    % best-fit model
```

The final values are reasonably close to the values $V=2$, $K_m=0.5$ that we used to generate the noise-free data. Graphically, the model looks close to the original data:

```{code-cell}
final_misfit_norm = norm(model(s) - w) 
hold on, fplot(model, [0, 6])
title('Michaelis-Menten fitting')    
labels = [labels, "nonlinear fit"];    
legend(labels, 'location', 'east');
```

For this particular model, we also have the option of linearizing the fit process. Rewrite the model as 

```{math}
:enumerated: false
\frac{1}{w} = \frac{\alpha}{s} + \beta = \alpha \cdot s^{-1} + \beta
```

for the new fitting parameters $\alpha=K_m/V$ and $\beta=1/V$. This corresponds to the misfit function whose entries are

$$f_i([\alpha,\beta]) = \left(\alpha \cdot \frac{1}{s_i} + \beta\right) - \frac{1}{w_i}$$

for $i=1,\ldots,m$. Although this misfit is nonlinear in $s$ and $w$, it's linear in the unknown parameters $\alpha$ and $\beta$. This lets us pose and solve it as a linear least-squares problem.

```{code-cell}
u = 1 ./ w;
A = [s.^(-1), s.^0];  
z = A \ u;
alpha = z(1);  beta = z(2);
```

The two fits are different because they do not optimize the same quantities.

```{code-cell}
linmodel = @(s) 1 ./ (beta + alpha ./ s);
final_misfit_linearized = norm(linmodel(s) - w)
fplot(linmodel, [0, 6])
labels = [labels, "linearized fit"];    
legend(labels, 'location', 'east');
```

The truly nonlinear fit is clearly better in this case. It optimizes a residual for the original measured quantity rather than a transformed one we picked for algorithmic convenience.
``````