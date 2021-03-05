# Fitting functions to data

```{index} interpolation; by polynomials
```

````{proof:example} Julia demo
:class: demo
{doc}`demos/fitting-tempinterp`
````

In {doc}`../linsys/polyinterp` we saw how a polynomial can be used to interpolate data—that is, derive a continuous function that evaluates to give a set of prescribed values. But interpolation may not be appropriate in many applications.

```{index} data fitting
```

In many cases we can get better results by relaxing the interpolation requirement. In the polynomial case this allows us to lower the degree of the polynomial, which limits the number of local max and min points. Let $(t_i,y_i)$ for $i=1,\ldots,m$ be the given points. We will represent the data by the polynomial

```{math}
:label: vanderls
y \approx f(t) = c_1 + c_2t + \cdots + c_{n-1} t^{n-2} + c_n t^{n-1},
```

with $n<m$. Just as in {eq}`vandersystem`, we can express a vector of $f$-values by a matrix-vector multiplication:

```{math}
:label: vandersystemrect
\begin{bmatrix}
f(t_1)                               \\
f(t_2)                               \\
f(t_3)                               \\
\vdots                               \\
f(t_m)
\end{bmatrix} =
\begin{bmatrix}
1      & t_1    & \cdots & t_1^{n-1} \\
1      & t_2    & \cdots & t_2^{n-1} \\
1      & t_3    & \cdots & t_3^{n-1} \\
\vdots & \vdots &        & \vdots    \\
1      & t_m    & \cdots & t_m^{n-1} \\
\end{bmatrix}
\begin{bmatrix}
c_1                                  \\
c_2                                  \\
\vdots                               \\
c_n
\end{bmatrix}.
```

```{index} Vandermonde matrix
```

````{proof:example} Julia demo
:class: demo
{doc}`demos/fitting-tempfit`
````

Note that $\mathbf{V}$ has the same structure as the Vandermonde matrix in {eq}`vandersystem` but is $m\times n$, thus taller than it is wide. It's impossible in general to satisfy $m$ conditions with $n<m$ variables, so rather than solving a linear system $\mathbf{V} \mathbf{c} = \mathbf{y}$, we have to find an approximation $\mathbf{V} \mathbf{c} \approx \mathbf{y}$. Below we specify precisely what is meant by this, but first we note that MATLAB conveniently uses the same backslash notation to solve the problem for both square and rectangular matrices.

## The least squares formulation

In the most general terms, our fitting functions take the form

```{math}
:label: fitbasis
f(t) = c_1 f_1(t) + \cdots + c_n f_n(t),
```

```{margin}
The essential feature of a linear least squares problem is that the fit depends only *linearly* on the unknown parameters.
```

where $f_1,\ldots,f_n$ are all known functions with no undetermined parameters. This leaves only $c_1,\ldots,c_n$ to be determined. The essential feature of a linear least squares problem is that the fit depends only *linearly* on the unknown parameters. For instance, a function of the form $f(t)=c_1 + c_2 e^{c_3 t}$ is not of this type.

At each observation $(t_i,y_i)$, we define a residual, $y_i - f(t_i)$. A sensible formulation of the fitting criterion is to minimize

```{math}
  R(c_1,\ldots,c_n) = \sum_{i=1}^m\, [ y_i - f(t_i) ]^2,
```

over all possible choices of parameters $c_1,\ldots,c_n$. We can apply linear algebra to write the problem in the form $R=\mathbf{r}^T \mathbf{r}$, where

```{math}
\mathbf{r} =
\begin{bmatrix}
y_1 \\ y_2 \\ \vdots \\y_{m-1} \\ y_m
\end{bmatrix} -
\begin{bmatrix}
f_1(t_1) & f_2(t_1) & \cdots & f_n(t_1) \\[1mm]
f_1(t_2) & f_2(t_2) & \cdots & f_n(t_2) \\[1mm]
& \vdots \\
f_1(t_{m-1}) & f_2(t_{m-1}) & \cdots & f_n(t_{m-1}) \\[1mm]
f_1(t_m) & f_2(t_m) & \cdots & f_n(t_m) \\[1mm]
\end{bmatrix}
\begin{bmatrix}
c_1 \\ c_2 \\ \vdots \\ c_n
\end{bmatrix}.
```

```{index} least squares
```

Recalling that $\mathbf{r}^T\mathbf{r}=\| \mathbf{r} \|_2^2$, and renaming the variables to standardize the statement, we arrive at the general {term}`linear least squares problem`:

````{proof:definition}
Given $\mathbf{A}\in\mathbb{R}^{m \times n}$ and $\mathbf{b}\in\mathbb{R}^m$, with $m>n$, find

```{math}
:label: linls
\operatorname{argmin}_{\mathbf{x}\in \mathbb{R}^n}  \| \mathbf{b}-\mathbf{A} \mathbf{x} \|_2^2.
```
````

The notation **argmin** above means to find an $\mathbf{x}$ that produces the minimum value.

```{index} norm; vector
```

While we could choose to minimize in any vector norm, the 2-norm is the most common and convenient choice. *For the rest of this chapter we exclusively use the 2-norm.* An $\mathbf{x}$ that achieves the minimum residual is considered to be a solution of the problem. In the edge case $m=n$ for a nonsingular $\mathbf{A}$, the definitions of the linear least squares and linear systems problems coincide: the solution of  $\mathbf{A}\mathbf{x}=\mathbf{b}$  implies $\mathbf{r}=\boldsymbol{0}$, which is a global minimum of $\| \mathbf{r} \|_2^2 \ge 0$.

## Change of variables

```{index} data fitting; by straight line
```

The most familiar and common case of a polynomial least squares fit is the straight line, $f(t) = c_1 + c_2 t$. Certain other fit functions can be transformed into this situation. For example, suppose we want to fit data using $g(t)= a_1 e^{a_2 t}$. Then

```{math}
:label: exptransfit
\log y \approx \log g(t) = (\log a_1) + a_2 t = c_1 + c_2 t.
```

```{index} data fitting; power law
```

While the fit of the $y_i$ to $f(t)$ is nonlinearly dependent on fitting parameters, the fit of $\log y_i$ to a straight line is a linear problem. Similarly, the power-law relationship $y\approx f(t)=a_1 t^{a_2}$ is equivalent to

```{math}
:label: powtransfit
\log y \approx (\log a_1) + a_2 (\log t).
```

````{proof:example} Julia demo
:class: demo
{doc}`demos/fitting-pirate`
````

Thus the variable $z=\log y$ can be fit linearly in terms of the variable $s=\log t$. In practice these two cases—exponential fit and power-law—are easily detected by using log-linear or log-log plots, respectively.

## Exercises

1. ✍ Suppose $f(x)$ is a differentiable nonnegative real function. Use calculus to prove that the local minimizers of $f(x)$ and $[f(x)]^2$ are the same values of $x$.

    (problem-fitcensus)=
2. ⌨  Here are counts of the U.S. population from the census performed every ten years, beginning with 1790 and ending with 2010.

    ``` julia
    3.929, 5.308, 7.240, 9.638, 12.87, 17.07, 23.19, 31.44, 39.82, 50.19, 62.95, 76.21,
    92.22, 106.0, 122.8, 132.2, 150.7, 179.3, 203.3, 226.5, 248.7, 281.4, 308.7
    ```

    **(a)** Find a best-fitting cubic polynomial for the data. Plot the data as points superimposed on a (smooth) graph of the cubic over the full range of time. Label the axes. What does the fit predict for the population in the years 2000,2010,2020?

    **(b)**Repeat (a) but use a fitting function of the form $y(t)\approx a e^{b t}$.

    **(c)**Look up the actual U.S. population in 2000 and 2010 and compare to the predictions of parts (a) and (b).

    % To make $y\approx a e^{bt}$ into a linear least squares problem, take the log (natural) of both sides: $\log(y) = \log(a) + b t$; then, let $c_1=\log(a)$ and $c_2=b$ and solve for the $c_i$.

    %% (a)
    <!-- load census;
    t = cdate - 1790;  y = pop;
    A = zeros(length(t),3);
    for j = 0:2
        A(:,3-j) = t.^j;
    end

    c = A\y;

    tt = linspace(0,t(end),200)';
    yy = polyval(c,tt);
    plot(t,y,'o',tt,yy)
    xlabel('year since 1790'), ylabel('population (millions)')

    pop2000 = polyval(c,210)
    pop2010 = polyval(c,220) -->

    %% (b)
    <!-- load census;
    t = cdate - 1790;  y = pop;
    A = [t.^1 t.^0];
    c = A\log(y);
    a = exp(c(2)), b = c(1)

    tt = linspace(0,t(end),200)';
    yy = a*exp(b*tt);
    plot(t,y,'o',tt,yy)
    xlabel('year since 1790'), ylabel('population (millions)')

    pop2000 = a*exp(210*b)
    pop2010 = a*exp(220*b) -->

3. ⌨  In this problem you are trying to find an approximation to the periodic function $f(t)=e^{\sin(t-1)}$ over one period, $0\le t \le 2\pi$. In Julia, let `t=LinRange(0,2pi,200)` and let `b` be a column vector of evaluations of $f$ at those values.

    **(a)** Find the coefficients of the least-squares fit
  
    ```{math}
      f(t) \approx c_1 + c_2t + \cdots + c_7 t^6.
    ```

     **(b)** Find the coefficients of the least-squares fit

    ```{math}
    f(t) \approx d_1 + d_2\cos(t) + d_3\sin(t) + d_4\cos(2t) +
    d_5\sin(2t).
    ```

     **(c)** Plot the original function $f(t)$ and the two approximations from (a) and (b) together on a well-labeled graph.
  
    %%  Problem 3.1.2
    % Try different fits on periodic function
    % y = exp( sin(t-1) )
    %% Part (a), polynomial fit
    %Use a polynomial fit to the data.
    <!-- clear all; format compact;
    t = linspace(0,2*pi,200)';  b = exp(sin(t-1));
    p = 7;  % degree of polynomial fit is p-1
    A = zeros(length(t),p+1);
    for j=1:p+1
        A(:,j) = t.^(j-1);
    end
    disp(' ')
    disp('Part (a) coefficients: ')
    c = (A'*A)\(A'*b)
    p = A*c;  % evaluate polynomial at all t values in one line! -->

    %% Part (b), trig function fit
    % Trigonometic function fit uses periodic functions which have that
    % property in common with the original data.
    <!-- A = ones(length(t),5);
    A(:,2) = cos(t);
    A(:,3) = sin(t);
    A(:,4) = cos(2*t);
    A(:,5) = sin(2*t);

    disp(' ')
    disp('Part (b) coefficients: ')
    d = (A'*A)\(A'*b)
    trig = A*d;  % again, compute the fits at all t values in one line! -->

    %% Part (c)
    % Plot the results for comparison.

    <!-- plot(t,p,'-',t,trig,'--',t,b,'-.')
    axis([0 2*pi 0 3])
    xlabel('t');
    legend('6th degree poly','trig fit','actual','Location','NorthEast');
    title('Problem 3.1.4 -- Fit of Periodic Function'); -->

    %% Discussion
    % The fit with the trigonometric functions is better near the ends
    % because those functions for the interpolation are periodic, like the
    % original function is.  The polynomial does pretty well, but to achieve a
    % specified level of error we expect that using a smaller number of trig
    % functions should do the job.

4. ⌨ Define the following data in Julia.

    ``` julia
    t = 0:.5:10
    y = tanh.(t)
    ```

     **(a)** Fit the data to a cubic polynomial. Plot the data together with the polynomial fit.

      **(b)** Fit the data to the function $c_1 + c_2z + c_3z^2 + c_4z^3$, where $z=t^2/(1+t^2)$. Plot the data together with the fit. What feature of $z$ makes this fit much better than the original cubic?
  
      <!-- t = (0:.5:10)';
      y = tanh(t);
      plot(t,y,'o')
      z = t.^2./(1+t.^2);
      A = [ t.^0 t t.^2 t.^3 ];
      c = A\y
      hold on
      plot(t,A*c) -->
      %% The key is that $z\to 1$ as $t\to \infty$.

5. ⌨  One series for finding $\pi$ is

    ```{math}
    \frac{\pi}{2} = 1 + \frac{1}{3} + \frac{1\cdot 2}{3\cdot5} + \frac{1\cdot 2\cdot 3}{3\cdot 5\cdot 7} + \cdots.
    ```

    Define $s_k$ to be the sum of the first $k$ terms (that is, up to $j=k=1$) and let $e_k=|\pi-2s_k|$. Do a numerical experiment to determine which is a better fit: $e_k\approx a b^k$, or $e_k\approx a k^b$. Then determine the parameters in that fit.

    <!-- a = 2;   % first term
    p = a;
    for j = 1:40
        a = a*j/(2*j+1);    % next term
        p(j+1) = p(j) + a;
    end
    e = abs(pi-p'); -->

    <!-- subplot(1,2,1),semilogy(e)
    subplot(1,2,2),loglog(e) -->

    %%
    % This suggests a fit $ab^k$.

    %%
    <!-- A = [ (0:40)', ones(41,1) ];
    c = A\log(e) -->

    %%
    <!-- b = exp(c(1))
    a = exp(c(2)) -->

6. ⌨  Kepler found that the orbital period $\tau$ of a planet depends on its mean distance $R$ from the sun according to $\tau=c R^{\alpha}$ for a simple rational number $\alpha$. Validate Kepler's result by using a linear least squares fit and the following table.

    | Planet  | Distance from sun in Mkm |  Orbital period in days              |
    |---------|:--------------------------------------------:|:-----------:|
    | Mercury | 57.59                                      | 87.99  |
    | Venus   | 108.11                                     | 224.7  |
    | Earth   | 149.57                                     | 365.26 |
    | Mars    | 227.84                                     | 686.98 |
    | Jupiter | 778.14                                     | 4332.4 |
    | Saturn  | 1427                                        | 10759   |
    | Uranus  | 2870.3                                     | 30684   |
    | Neptune | 4499.9                                     | 60188   |
    | Pluto   | 5909                                        | 90710   |

7. ✍ Show that finding a best fit of the form

    ```{math}
    y(t) \approx \frac{a}{t+b}
    ```

    can be transformed into a linear fitting problem (with different undetermined coefficients) by rewriting the equation.

    % We have $1/y \approx (t/a) + (a/b) = ct +d$, which is a straight line fitting problem.

8. ✍ Show how to find the constants $a$ and $b$ in a data fitting problem of the form $y(t)\approx t/(at+b)$ by transforming it into a linear least squares fitting problem.
