# Nonlinear least squares

After the solution of square linear systems, we generalized to the case of having more constraints to satisfy than available variables. Our next step is to do the same for nonlinear equations, thus filling out this table:

|                    | **linear** | **nonlinear** |
|:---------------:|:---------:|:--------------:|
| **square**         | $\mathbf{A}\mathbf{x}=\mathbf{b}$ |$\mathbf{f}(\mathbf{x})=\boldsymbol{0}$ |
| **overdetermined** | $\min\, \bigl\|\mathbf{A}\mathbf{x} - \mathbf{b}\bigr\|_2$ | $\min\, \bigl\|\mathbf{f}(\mathbf{x}) \bigr\|_2$ |

```{index} nonlinear least squares
```

We now define the {term}`nonlinear least squares problem`.

```{prf:definition} Nonlinear least squares problem
Given a function $\mathbf{f}(\mathbf{x})$ mapping from $\real^n$ to $\real^m$, find $\mathbf{x}\in\real^n$ such that $\bigl\|\mathbf{f}(\mathbf{x})\bigr\|_2$ is minimized.
```

As in the linear case, we consider only overdetermined problems, where $m>n$. Minimizing a positive quantity is equivalent to minimizing its square, so we could also define the result as minimizing $\phi(\mathbf{x})=\mathbf{f}(\mathbf{x})^T\mathbf{f}(\mathbf{x})$.

```{index} Gauss--Newton method
```

You should not be surprised to learn that we can formulate an algorithm by substituting a linear model function for $\mathbf{f}$. At a current estimate $\mathbf{x}_k$ we define

```{math}
  \mathbf{q}(\mathbf{x})  = \mathbf{f}(\mathbf{x}_k) + \mathbf{A}_k(\mathbf{x}-\mathbf{x}_k),
```

where $\mathbf{A}_k$ might be the exact $m\times n$ Jacobian matrix, $\mathbf{J}(\mathbf{x}_k)$, or a finite-difference or Broyden approximation of it as described in {doc}`quasinewton`.

In the square case we solved $\mathbf{q}=\boldsymbol{0}$ to define the new value for $\mathbf{x}$, leading to the condition $\mathbf{A}_k\mathbf{s}_k=-\mathbf{f}_k$, where  $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$. Now, with $m>n$, we cannot expect to solve $\mathbf{q}=\boldsymbol{0}$, so instead we define $\mathbf{x}_{k+1}$ as the value that minimizes $\| \mathbf{q} \|_2$:

```{math}
  :label: NLSmin
  \mathbf{x}_{k+1} = \mathbf{x}_k + \mathbf{s}_k, \text{ where } \| \mathbf{A}_k \mathbf{s}_k + \mathbf{f}(\mathbf{x}_k)\|_2 \text{ is minimized.}
```

```{margin}
Gauss--Newton solves a series of linear least-squares problems in order to solve a nonlinear least squares problem.
```

We have just described the {term}`Gauss--Newton method`. Gauss--Newton solves a series of linear least-squares problems in order to solve a nonlinear least squares problem.

Here we reap an amazing benefit from Julia: The functions {ref}`function-newtonsys` and {ref}`function-levenberg`, which were introduced for the case of $m=n$ nonlinear equations, work *without modification* as the Gauss--Newton method for the overdetermined case! The reason is that the backslash operator applies equally well to the linear system and linear least squares problems, and nothing else in the function was written with explicit reference to $n$.

## Convergence

```{margin}
Gauss--Newton converges less rapidly than Newton's method if the optimal residual is not small.
```

In the multidimensional Newton method for a nonlinear system, we expect quadratic convergence to a solution in the typical case. For the Gauss--Newton method, the picture is more complicated. As always in least squares problems, the residual $\mathbf{f}(\mathbf{x})$ will not necessarily be zero when $\|\mathbf{f}\|$ is minimized. Suppose that the minimum value of $\|\mathbf{f}\|$ is $R>0$. In general, we might observe quadratic-like convergence until the iterate $\|\mathbf{x}_k\|$ is within distance $R$ of a true minimizer, and linear convergence thereafter. When $R$ is not sufficiently small, the convergence can be quite slow.

## Nonlinear data fitting

```{index} data fitting; nonlinear
```

In {doc}`../leastsq/fitting` we saw how to fit functions to data values, provided that the set of candidate fitting functions depends linearly on the undetermined coefficients. We now have the tools to generalize that process to fitting functions that depend nonlinearly on unknown parameters. Suppose that $(t_i,y_i)$ for $i=1,\ldots,m$ are given points. We wish to model the data by a function $g(t;\mathbf{c})$ that depends on unknown parameters $c_1,\ldots,c_n$ in an arbitrary way. A standard approach is to minimize the discrepancy between the model and the observations, in a least squares sense:

```{math}
  \min_{\mathbf{c}\in\mathbb{R}^n} \sum_{i=1}^m \bigl[ g(t_i;\mathbf{c})-y_i \bigr]^2 = \min_{\mathbf{c}\in\mathbb{R}^n} \bigl\| \mathbf{f}(\mathbf{c}) \bigr\|^2,
```

where $\mathbf{f}(\mathbf{c})$ is the vector of values $g(t_i;\mathbf{c})-y_i$. We call $\mathbf{f}$ a **misfit** function: the smaller the norm of the misfit, the better the fit.

````{prf:example} Julia demo
:class: demo
:label: demos-nlsq-MM
{doc}`demos/nlsq-MM`
````

The form of $g$ is up to the modeler. There may be compelling theoretical choices, or you may just be looking for enough algebraic power to express the data well. Naturally, in the special case where the dependence on $\mathbf{c}$ is linear, i.e.,

```{math}
  g(t;\mathbf{c}) = c_1 g_1(t) + c_2 g_2(t) + \cdots + c_m g_m(t),
```

then the misfit function is also linear in $\mathbf{c}$ and the fitting problem reduces to linear least squares.

## Exercises

(problem-NLSonevar)=

1. ✍ Define $\mathbf{f}(x)=[ x-8, \; x^2-4 ]^T$.

    **(a)** Write out the linear model of $\mathbf{f}$ at $x=2$.

    **(b)** Find the estimate produced by one step of the Gauss--Newton method, starting at $x=2$.
  
    ````{only} solutions

    %% Problem 4.5.1
    %
    % (a) We want to find an approximation to the $$ x $$ that minimizes the
    % two norm of $$ f(x) $$.  The linearization is analogous to Newton's
    % method for doing this.  We have that $$ f(x) = [ x-8; x^2-4 ] $$.
    %
    % The linearization is
    %
    % $$ L(x) = f(2) - \frac{df}{dx}(2) (x-2) $$.
    %
    % The Jacobian is $$ df/dx =  [1; 2x] $$, so that the the linearization is
    %
    % $$ L(x) = \left[ \begin{array}{c} -6 \\ 0 \end{array}  \right] + \left[ \begin{array}{c} 1 \\ 4 \end{array}  \right] (x - 2) $$.
    %
    % (b) Taking one step means solving the normal equations to minimize the
    % 2-norm of $$ L(x)=0 $$.  The normal equations are
    % found by premultiplying by the transpose of the Jacobian, $$ [1 4] $$,
    % and then one solves for $$ x $$. One obtains
    %
    % $$ -[1\ 4]\left[ \begin{array}{c} -6 \\ 0 \end{array}  \right] = [1\ 4]\left[ \begin{array}{c} 1 \\ 4 \end{array}  \right] (x-2) \Rightarrow 6 = 17(x-2) $$
    %
    % One could solve for $$ x-2 = 6/17 $$ and add the result to 2, or solve
    % directly for $$ x = 2+6/17 = 40/17 $$.

    ````

2. ✍ (continuation of preceding problem) The Gauss--Newton method replaces $\mathbf{f}(\mathbf{x})$ by a linear model and minimizes the residual norm of it. An alternative is to replace $\| \mathbf{f}(\mathbf{x}) \|_2^2$ by a scalar *quadratic* model $q(\mathbf{x})$, and minimize that.
  
    **(a)** Using $\mathbf{f}(x) = [ x-8, \; x^2-4 ]^T$, let $q(x)$ be defined by the first three terms in the Taylor series for $\| \mathbf{f}(x) \|_2^2$ at $x=2$.

    **(b)** Find the unique $x$ that minimizes $q(x)$. Is the result the same as the estimate produced by Gauss--Newton?
  
3. ⌨  A famous result by Kermack and McKendrick in 1927 {cite}`kermackContributionMathematical1927` suggests that in epidemics that kill only a small fraction of a susceptible population, the death rate as a function of time is well modeled by

    ```{math}
    w'(t) = A \operatorname{sech}^2[B(t-C)],
    ```
  
    for constant values of the parameters $A,B,C$. Since the maximum of sech is $\operatorname{sech}(0)=1$, $A$ is the maximum death rate and $C$ is the time of peak deaths. You will use this model to fit the deaths per week from plague recorded in Mumbai in a period during 1906:

    ``` julia
    5, 10, 17, 22, 30, 50, 51, 90, 120, 180, 292, 395, 445, 775, 780,
    700, 698, 880, 925, 800, 578, 400, 350, 202, 105, 65, 55, 40, 30, 20
    ```

    **(a)** Use {ref}`function-levenberg` to find the best least-squares fit to the data using the $\operatorname{sech}^2$ model. Make a plot of the model fit superimposed on the data. What are $A$ and $C$?

    **(b)** In practice, one would like a model to predict the full course of the epidemic before it has reached its peak and subsided. Redo the fitting from part~(a) using only the first 12 data values. Add this model to your plot and report $A$ and $C$. Is this model a useful predictor of the value and timing of the maximum death rate?

    **(c)** Repeat part (b) using the first 13 data values.

    ````{only} solutions
    ``` matlab
    %%
    % A famous result by Kermack and McKendrick in 1927 suggests that in
    % epidemics that kill only a small fraction of a susceptible population,
    % the number of deaths as a function of time is well modeled by
    %
    %  $$w'(t) = A \sech^2(B(t-C))$$
    %
    % for constant values of the parameters $A,B,C$. We will use this to fit
    % the deaths per week from plague recorded in Bombay (now Mumbai) in a
    % period during 1906.

    deaths=[ ...
        5, 10, 17, 22, 30, 50, 51, 90, 120, 180, 292, 395, 445, 775, 780, ...
        700, 698, 880, 925, 800, 578, 400, 350, 202, 105, 65, 55, 40, 30, 20 ]';

    %%
    % Here is the residual function for the fitting.
    m = length(deaths);
    week = (1:m)';
    model = @(x,t) x(1)*sech(x(2)*(t-x(3))).^2;
    f = @(x) model(x,week(1:12)) - deaths(1:12);

    %%
    x = levenberg(f,[1000;1;15]);
    ABC = x(:,end);
    A = ABC(1),  B = ABC(2),  C = ABC(3)

    %%
    % The model projects a maximum of about 880 deaths per week.

    %%
    hold on
    fplot( @(t) model(x(:,end),t), [0 30])
    ```
    ````

4. ⌨  The *Rosenbrock banana function* is defined as $\| \mathbf{f}(\mathbf{x}) \|_2^2$, where
  
    ```{math}
    \mathbf{f}(\mathbf{x}) =
    \begin{bmatrix}
      10(x_2-x_1^2) \\ 1-x_1
    \end{bmatrix}.
    ```

    Use {ref}`function-newtonsys` to find a minimizer of the banana function starting from $(-1.4,5.1)$. (If you're curious about its name, make a contour plot of the residual over $-2\le x_1 \le 3$, $-1\le x_2 \le 4$.) Show all the digits of the final result.

5. ⌨  In this problem you are to fit a function of the form
  
    ```{math}
    P(t) = c_1 + c_2 e^{c_3 t}
    ```

    to US census data for the twentieth century. Starting in 1900, the population in millions every ten years was:
  
    ``` julia
    76.0, 92.0, 105.7, 122.8, 131.7, 150.7, 179.0, 205.0, 226.5, 248.7
    ```
  
    Use nonlinear least squares to determine the unknown parameters $c_1$, $c_2$, $c_3$ in $P$. To aid convergence, rescale the data using the time variable $t = (\text{year}-1900)/100$ and divide the population numbers above by $100$. Using your model $P(t)$, predict the result of the 2000 census, and compare it to the exact figure (which can be found easily on the internet).

    ````{only} solutions

    ``` matlab
    %% Problem 4.5.3
    %
    % We need to enter and scale the data; that happens first.  Then create
    % scalar function f and its Jacobian (which is rectangular here).
    % Then, send them to newtonsys and plot the result.

    % create data
    year = [1900, 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990]';
    pop = [76.0, 92.0, 105.7, 122.8, 131.7, 150.7, 179.0, 205.0, 226.5, 248.7]';
    t = (year-year(1))/100;   % better independent variable
    y = pop/100;

    % function to minimize and Jacobian
    f = @(x) y-x(1)-x(2)*exp(x(3)*t);
    Jac = @(x) [-ones(size(y)), -exp(x(3)*t), (-x(2)*t).*exp(x(3)*t)];

    % Now use newtonsys.m
    x0 = [7e-1; 1e1; 0.1];  % initial guess
    x = newtonsys(f,Jac,x0);
    c = x(:,end)

    popest = c(1)+c(2)*exp(c(3)*t);              % the fit
    plot(year,pop,'bo',year,popest*100,'-b','LineWidth',2)  % plot it
    xlabel('Year'), ylabel('Population (in millions)')
    legend('Data','Fit','Location','NorthWest')
    legend('boxoff')

    % 2000 population estimate in millions follows; this would be at t=1.
    format short
    pop2000est = 1e2*(c(1)+c(2)*exp(c(3)*1))

    % Population estimate for 2000 from web is 282.2 million.
    ```
    ````

6. ⌨ The position of the upper lid during an eye blink can be measured from high speed video {cite}`wuEffectsMild2014`, and it may be possible to classify blinks based in part on fits to the lid position {cite}`broschBlinkCharacterization2017`. The lid position functions proposed to fit blinks is a product of a monomial or polynomial multiplying a decaying exponential {cite}`berkeKineticsLid1998`.  In this problem, you will generate representative data, add a small amount of noise to it, and then perform nonlinear least squares fits to the data.

   **(a)** Consider the function $y(\mathbf{a}) = a_1 \exp \left( -a_2 t^{a_3} \right)$, using the vector of coefficients $\mathbf{a} = [a_1 \ a_2\ a_3]$, and create eyelid position data as follows:

    ``` matlab
    N = 20;                            % number of time values
    t = linspace(1,N,N)'/N;            % equally spaced to t=1
    a = [10, 10, 2]';                  % baseline values
    y = a(1)*t.^2.*exp(-a(2)*t.^a(3)); % ideal data
    ym = y;                            % vector for data
    ir = 1:N-1;                        % range to add noise
    rng(13);                           % set seed for rand
    noise = 0.03;                      % amplitude of noise
    ym(ir) = y(ir) + noise*rand(size(y(ir)));  % add noise
    ```

    **(b)** Using the data `(t,ym)`, find the nonlinear least squares fit using {ref}`function-levenberg`.

    **(c)** Plot the fits using `np = 100` points over `t=(1:np)/np` together with symbols for the `N` measured data points `ym`.

    **(d)** Increase the noise to 5% and 10%. You may have to increase the number of measured points `N` and/or the maximum number of iterations.  How close are the coefficients?  Plot the data and the resulting fit for each case.

    ````{only} solutio

    function NLLS_blink_fit_ex1_fun
    %% model blink fit problem
    % Based fitting from Brosch et al, J Modeling Ophthalmol 2017
    % using data from Ziwei Wu et al, IOVS, 2015
    %
    %% setup
    close all, clear all, clc
    format short e, format compact
    rng(13)
    %% create synthetic blink data
    %
    N = 20; noise = 0.1;
    t = [1:N]'/N;
    a = [20, -10, -8, 7, 2]';  % reasonable initial guess
    y = (a(1)+a(2)*t+a(3)*t.^2).*t.^2.*exp(-a(4)*t.^a(5));
    plot(t,y,'LineWidth',2)
    xlabel('t'), ylabel('y')
    ym = y;
    irange = 1:N;
    ym(irange) = y(irange) + noise*rand(size(y(irange)));
    hold on
    plot(t,ym,'or','LineWidth',2)
    legend('base fn','data','Location','Best')
    hold off

    %% Using Matlab builtin
    %
    % Use fminsearch to solve the nonlinear set of equations
    % This method does not need a Jacobian (uses Nelder-Mead)
    %

    % Sum of squares function
    fblink = @(x) sum(( ym-(x(1)+x(2)*t+x(3)*t.^2).*t.^2.*exp(-x(4)*t.^x(5)) ).^2);
    x0 = a;
    myopts = optimset('MaxFunEvals',2000,'MaxIter',2000)
    [x2,fval2] = fminsearch(@(x) fblink(x),x0,myopts)
    Nplot = 100;
    tt = linspace(1,Nplot,Nplot)'/Nplot;
    yhat = (x2(1)+x2(2)*tt+x2(3)*tt.^2).*tt.^2.*exp(-x2(4)*tt.^x2(5));
    figure
    plot(t,ym,'ro',tt,yhat,'--g','LineWidth',2)
    legend('Data','fminsearch','Location','Best');

    figure
    yy = (a(1)+a(2)*tt+a(3)*tt.^2).*tt.^2.*exp(-a(4)*tt.^a(5));
    plot(tt,yy,'-b',tt,yhat,'--g','LineWidth',2)

    %% Levenberg-Marquardt type approach
    %
    % Uses LevMarqRB1.m of text with default tolerances and max iterations
    %

    lambda = 10;   % iteration parameter
    % Using the following function allows passing of t and y as parameters
    % Equations formulated like normal equations in this function
    LMmaxit = 500;    % max allowed iterations
    [xans,fans,err,iter] = LevMarqRB1(@fblinkLM,x0,lambda,1e4*eps,LMmaxit,t,ym);
    if iter>=LMmaxit
        iter
        disp('Max iterations exceeded for L-M method')
    end
    x3 = xans(:,end)
    yhat3 = (x3(1)+x3(2)*tt+x3(3)*tt.^2).*tt.^2.*exp(-x3(4)*tt.^x3(5));

    figure
    plot(t,ym,'ro',tt,yhat3,'--g','LineWidth',2)
    legend('Data','Lev-Marq','Location','Best');

    figure
    semilogy(err,'LineWidth',2)
    xlabel('iteration number'), ylabel('residual')

    function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % LevMarqRB1: Levenberg-Marquardt-like algorithm with adaptive lambda
    %    [x,f,ea,iter]=newtmult(func,x0,es,maxit,p1,p2,...):
    %     uses the Newton-Raphson method to find the roots of
    %     a system of nonlinear equations
    %    Modified version from Chapra to accumulate iterates, function values,
    %    and approximate errors; for exploration of method
    % input:
    %     func = name of function that returns f and J
    %     x0 = initial guess (column vector)
    %     lambda = initial iteration parameter (scalar; default = 25 if absent)
    %     errtol = desired percent relative error (scalar; default = 0.0001%)
    %     maxit = maximum allowable iterations (default = 50)
    %     p1,p2,... = additional parameters used by function
    % output:
    %     x = vector of iterates
    %     f = vector of functions evaluated at iterates; returns list at end
    %     ea = approximate percent relative error (%) (really % change in iterates)
    %     iter = number of iterations

    if nargin<2,error('at least 2 input arguments required'),end
    if nargin<3|isempty(lambda),lambda=25;end
    if nargin<4|isempty(errtol),errtol=1e-6;end
    if nargin<5|isempty(maxit),maxit=50;end
    iter = 0;
    x=x0;
    fiter = [];
    xiter = x0;
    ea = [100];
    while (ea(end) > errtol & iter < maxit)
      [f,J]=func(x,varargin{:});         %  name supplied as input
      JtJ = J'*J;
      A = JtJ+lambda*diag(diag(JtJ));    %  matrix for L-M method
      b = -(J'*f);                       %  rhs for L-M
      dx = A\b;                          % solve for update
      x = x+dx;
      fiter = [fiter f];
      xiter = [xiter x];
      iter = iter + 1;
      ea=[ea; 100*norm(dx./x,inf)]; % relative percent change in iterate
      % add a naive adjustment to L-M of text
      if ea(end)-ea(end-1)<0        % probably different for each problem
          lambda = lambda/10;       % reduce lambda if iterates get very close
      else
          lambda = lambda*2;
      end
    end
    [f,J]=func(x,varargin{:});
    f = [fiter f];
    x = xiter;
    end

    function [f,Jac] = fblinkLM(x,t,y)
    % function [f,Jac] = fblinkLM(x,t,y)
    % input:
    %    x = independent variable, column vector of length 5
    %    a = parameter for ellipse semi axis, scalar
    % output:
    %   f = function values, column vector of length m
    %   Jac = Jacobian matrix mx2

    f(:,1) = y-(x(1)+x(2)*t+x(3)*t.^2).*t.^2.*exp(-x(4)*t.^x(5));

    Jac(:,1) = -t.^2.*exp(-x(4)*t.^x(5));
    Jac(:,2) = -t.^3.*exp(-x(4)*t.^x(5));
    Jac(:,3) = -t.^4.*exp(-x(4)*t.^x(5));
    Jac(:,4) = (x(1)+x(2)*t+x(3)*t.^2).*t.^2.*exp(-x(4)*t.^x(5)).*t.^x(5);
    Jac(:,5) = (x(1)+x(2)*t+x(3)*t.^2).*t.^2.*exp(-x(4)*t.^x(5))*...
                x(4).*t.^x(5).*log(t);
    end
    end
    ````

7. ⌨ Repeat the previous problem using the fitting function $y(\mathbf{a}) = (a_1+a_2 t + a_3 t^2) t^2 \exp \left( -a_4 t^{a_5} \right)$, using the vector of coefficients $\mathbf{a} = [a_1 \ a_2\ a_3\ a_4\ a_5]$. (This was the choice used in Brosch et al {cite}`broschBlinkCharacterization2017`.)  Use `a = [20, -10, -8, 7, 2]` to create the data and as an initial guess for the coefficients for the fit to the noisy data.

    ````{only} solutions
    
    function NLLS_blink_fit_ex2_fun
    %% model blink fit problem
    % Based fitting from Brosch et al, J Modeling Ophthalmol 2017
    % using data from Ziwei Wu et al, IOVS, 2015
    %
    %% setup
    close all, clear all, clc
    format short e, format compact
    rng(13)
    %% create synthetic blink data
    %
    N = 20; noise = 0.03;
    t = linspace(1,N,N)'/N;
    a = [20, 10, 2]';  % reasonable initial guess
    y = a(1)*t.^2.*exp(-a(2)*t.^a(3));
    plot(t,y,'LineWidth',2)
    xlabel('t'), ylabel('y')
    ym = y;
    irange = 1:N-1;
    ym(irange) = y(irange) + noise*rand(size(y(irange)));
    hold on
    plot(t,ym,'or','LineWidth',2)
    legend('base fn','data','Location','Best')
    hold off
    %
    % use initial guesses of b,then use nonlinear least squares to find
    % new coefficients

    %% Using Matlab builtin
    %
    % Use fminsearch to solve the nonlinear set of equations
    % This method does not need a Jacobian (uses Nelder-Mead)
    %

    fblink = @(x) sum(( ym-x(1)*t.^2.*exp(-x(2)*t.^x(3)) ).^2);  % Note the sum of squares
    x0 = a;
    myopts = optimset('MaxFunEvals',2000,'MaxIter',1000)
    [x2,fval2] = fminsearch(@(x) fblink(x),x0,myopts)
    Nplot = 100;
    tt = linspace(1,Nplot,Nplot)'/Nplot;
    yhat = x2(1)*tt.^2.*exp(-x2(2)*tt.^x2(3));
    figure
    plot(t,ym,'ro',tt,yhat,'--g','LineWidth',2)
    legend('Data','fminsearch','Location','Best');

    figure
    yy = a(1)*tt.^2.*exp(-a(2)*tt.^a(3));
    plot(tt,yy,'-b',tt,yhat,'--g','LineWidth',2)

    %% Levenberg-Marquardt type approach
    %
    % Uses LevMarqRB1.m of text with default tolerances and max iterations
    %

    lambda = 10;   % iteration parameter
    % Using the following function allows passing of t and y as parameters
    % Equations formulated like normal equations in this function
    LMmaxit = 300;    % max allowed iterations
    [xans,fans,err,iter] = LevMarqRB1(@fblinkLM2,x0,lambda,1e4*eps,LMmaxit,t,ym);
    if iter>=LMmaxit
        iter
        disp('Max iterations exceeded for L-M method')
    end
    x3 = xans(:,end)
    yhat3 = x3(1)*tt.^2.*exp(-x3(2)*tt.^x3(3));

    figure
    plot(t,ym,'ro',tt,yhat3,'--g','LineWidth',2)
    legend('Data','Lev-Marq','Location','Best');

    figure
    semilogy(err,'LineWidth',2)
    xlabel('iteration number'), ylabel('residual')

    function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % LevMarqRB1: Levenberg-Marquardt-like algorithm with adaptive lambda
    %    [x,f,ea,iter]=newtmult(func,x0,es,maxit,p1,p2,...):
    %     uses the Newton-Raphson method to find the roots of
    %     a system of nonlinear equations
    %    Modified version from Chapra to accumulate iterates, function values,
    %    and approximate errors; for exploration of method
    % input:
    %     func = name of function that returns f and J
    %     x0 = initial guess (column vector)
    %     lambda = initial iteration parameter (scalar; default = 25 if absent)
    %     errtol = desired percent relative error (scalar; default = 0.0001%)
    %     maxit = maximum allowable iterations (default = 50)
    %     p1,p2,... = additional parameters used by function
    % output:
    %     x = vector of iterates
    %     f = vector of functions evaluated at iterates; returns list at end
    %     ea = approximate percent relative error (%) (really % change in iterates)
    %     iter = number of iterations

    if nargin<2,error('at least 2 input arguments required'),end
    if nargin<3|isempty(lambda),lambda=25;end
    if nargin<4|isempty(errtol),errtol=1e-6;end
    if nargin<5|isempty(maxit),maxit=50;end
    iter = 0;
    x=x0;
    fiter = [];
    xiter = x0;
    ea = [100];
    while (ea(end) > errtol & iter < maxit)
      [f,J]=func(x,varargin{:});         %  name supplied as input
      JtJ = J'*J;
      A = JtJ+lambda*diag(diag(JtJ));    %  matrix for L-M method
      b = -(J'*f);                       %  rhs for L-M
      dx = A\b;                          % solve for update
      x = x+dx;
      fiter = [fiter f];
      xiter = [xiter x];
      iter = iter + 1;
      ea=[ea; 100*norm(dx./x,inf)]; % relative percent change in iterate
      % add a naive adjustment to L-M of text
      if ea(end)-ea(end-1)<0        % probably different for each problem
          lambda = lambda/10;       % reduce lambda if iterates get very close
      else
          lambda = lambda*2;
      end
    end
    [f,J]=func(x,varargin{:});
    f = [fiter f];
    x = xiter;

    end
    function [f,Jac] = fblinkLM2(x,t,y)
    % function [f,Jac] = fblinkLM(x,t,y)
    % input:
    %    x = independent variable, column vector of length 5
    %    a = parameter for ellipse semi axis, scalar
    % output:
    %   f = function values, column vector of length m
    %   Jac = Jacobian matrix mx2

    f(:,1) = y-x(1)*t.^2.*exp(-x(2)*t.^x(3));

    Jac(:,1) = -t.^2.*exp(-x(2)*t.^x(3));
    Jac(:,2) = x(1)*t.^2.*exp(-x(2)*t.^x(3)).*t.^x(3);
    Jac(:,3) = x(1)*t.^2.*exp(-x(2)*t.^x(3)).*...
                x(2).*t.^x(3).*log(t);
    end
    end
    ````


    <!-- \item ⌨ The following problems ask you to generate nonlinear fits using the specified functions with the same levels of noise and numbers of points as in the previous two problems; however, use "randn" rather than "rand". For these functions, estimate to two digits the amplitude of the noise where the SSE becomes larger than unity.  For all computations, use initial guess $\mathbf{a} = [ 2,\ 2\pi, \  2.5]$.

   	1. 
    \item $y(\mathbf{a}) = a_1 t \sin \left( a_2 t^{a_3} \right)$.
    \item $y(\mathbf{a}) = a_1 (1-t) \cos \left( a_2 t^{a_3} -->

    ````{only} solutions
    function NLLS_blink_fit_ex3_fun
    %% model blink fit problem
    % Based fitting from Brosch et al, J Modeling Ophthalmol 2017
    % using data from Ziwei Wu et al, IOVS, 2015
    %
    %% setup
    close all, clear all, clc
    format short e, format compact
    rng(13)
    %% create synthetic blink data
    %
    N = 20; noise = 0.21;
    t = linspace(1,N,N)'/N;
    a = [2, 2*pi, 2.5]';  % reasonable initial guess
    y = a(1)*t.*sin( a(2)*t.^a(3) );
    plot(t,y,'LineWidth',2)
    xlabel('t'), ylabel('y')
    ym = y;
    irange = 1:N-1;
    ym(irange) = y(irange) + noise*randn(size(y(irange)));
    hold on
    plot(t,ym,'or','LineWidth',2)
    legend('base fn','data','Location','Best')
    hold off
    %
    % use initial guesses of b,then use nonlinear least squares to find
    % new coefficients

    %% Using Matlab builtin
    %
    % Use fminsearch to solve the nonlinear set of equations
    % This method does not need a Jacobian (uses Nelder-Mead)
    %

    fblink = @(x) sum(( ym-x(1)*t.*sin(x(2)*t.^x(3)) ).^2);  % Note the sum of squares
    x0 = a;
    myopts = optimset('MaxFunEvals',2000,'MaxIter',1000)
    [x2,fval2] = fminsearch(@(x) fblink(x),x0,myopts)
    Nplot = 100;
    tt = linspace(1,Nplot,Nplot)'/Nplot;
    yhat = x2(1)*tt.*sin(x2(2)*tt.^x2(3));
    figure
    plot(t,ym,'ro',tt,yhat,'--g','LineWidth',2)
    legend('Data','fminsearch','Location','Best');

    figure
    yy = a(1)*tt.*sin( a(2)*tt.^a(3) );
    plot(tt,yy,'-b',tt,yhat,'--g','LineWidth',2)

    %% Levenberg-Marquardt type approach
    %
    % Uses LevMarqRB1.m of text with default tolerances and max iterations
    %

    lambda = 10;   % iteration parameter
    % Using the following function allows passing of t and y as parameters
    % Equations formulated like normal equations in this function
    LMmaxit = 300;    % max allowed iterations
    [xans,fans,err,iter] = LevMarqRB1(@fblinkLM3,x0,lambda,1e4*eps,LMmaxit,t,ym);
    if iter>=LMmaxit
        iter
        disp('Max iterations exceeded for L-M method')
    end
    x3 = xans(:,end)
    fblink(x3)
    yhat3 = x3(1)*tt.*sin(x3(2)*tt.^x3(3));

    figure
    plot(t,ym,'ro',tt,yhat3,'--g','LineWidth',2)
    legend('Data','Lev-Marq','Location','Best');

    figure
    semilogy(err,'LineWidth',2)
    xlabel('iteration number'), ylabel('residual')

    function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % function [x,f,ea,iter]=LevMarqRB1(func,x0,lambda,errtol,maxit,varargin)
    % LevMarqRB1: Levenberg-Marquardt-like algorithm with adaptive lambda
    %    [x,f,ea,iter]=newtmult(func,x0,es,maxit,p1,p2,...):
    %     uses the Newton-Raphson method to find the roots of
    %     a system of nonlinear equations
    %    Modified version from Chapra to accumulate iterates, function values,
    %    and approximate errors; for exploration of method
    % input:
    %     func = name of function that returns f and J
    %     x0 = initial guess (column vector)
    %     lambda = initial iteration parameter (scalar; default = 25 if absent)
    %     errtol = desired percent relative error (scalar; default = 0.0001%)
    %     maxit = maximum allowable iterations (default = 50)
    %     p1,p2,... = additional parameters used by function
    % output:
    %     x = vector of iterates
    %     f = vector of functions evaluated at iterates; returns list at end
    %     ea = approximate percent relative error (%) (really % change in iterates)
    %     iter = number of iterations

    if nargin<2,error('at least 2 input arguments required'),end
    if nargin<3|isempty(lambda),lambda=25;end
    if nargin<4|isempty(errtol),errtol=1e-6;end
    if nargin<5|isempty(maxit),maxit=50;end
    iter = 0;
    x=x0;
    fiter = [];
    xiter = x0;
    ea = [100];
    while (ea(end) > errtol & iter < maxit)
      [f,J]=func(x,varargin{:});         %  name supplied as input
      JtJ = J'*J;
      A = JtJ+lambda*diag(diag(JtJ));    %  matrix for L-M method
      b = -(J'*f);                       %  rhs for L-M
      dx = A\b;                          % solve for update
      x = x+dx;
      fiter = [fiter f];
      xiter = [xiter x];
      iter = iter + 1;
      ea=[ea; 100*norm(dx./x,inf)]; % relative percent change in iterate
      % add a naive adjustment to L-M of text
      if ea(end)-ea(end-1)<0        % probably different for each problem
          lambda = lambda/10;       % reduce lambda if iterates get very close
      else
          lambda = lambda*2;
      end
    end
    [f,J]=func(x,varargin{:});
    f = [fiter f];
    x = xiter;
    end

    function [f,Jac] = fblinkLM3(x,t,y)
    % function [f,Jac] = fblinkLM(x,t,y)
    % input:
    %    x = independent variable, column vector of length 5
    %    a = parameter for ellipse semi axis, scalar
    % output:
    %   f = function values, column vector of length m
    %   Jac = Jacobian matrix mx2

    f(:,1) = y-x(1)*t.*sin(x(2)*t.^x(3));

    Jac(:,1) = -t.*sin(x(2)*t.^x(3));
    Jac(:,2) = -x(1)*t.*cos(x(2)*t.^x(3)).*t.^x(3);
    Jac(:,3) = -x(1)*t.*cos(x(2)*t.^x(3)).*...
                x(2).*t.^x(3).*log(t);
    end
    end
    ````
