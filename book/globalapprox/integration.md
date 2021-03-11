# Spectrally accurate integration
\label{sec:specquad}

In \charef{localfuncapprox} we derived methods of order 2, 4, and higher for numerical integration. (Recall that *quadrature* is another common term for numerical integration.) These all started with an expression of the form

$$
  :label: globalquad
  \int_{-1}^1 f(x)\, dx \approx \sum_{k=0}^n c_k f(t_k),
$$

for a collection of nodes $t_0,\ldots,t_n$ in $[-1,1]$ and weights $c_0,\ldots,c_n$. (Throughout this section we use $[-1,1]$ as the domain of the integral; for a general interval $[a,b]$, see {ref}`prob-specquad-integrateinterval`.) The nodes and weights are independent of the integrand $f(x)$ and determine the implementation and properties of the formula.

The process for deriving a specific method was to interpolate the integrand, then integrate the interpolant. Piecewise linear interpolation at equally spaced nodes, for instance, leads to the trapezoid formula. \texthighlight{specint}{When the integrand is approximated by a spectrally accurate global function, the integration formulas are also 
```{index} convergence rate!spectral
```
 spectrally accurate.}



## Periodic functions

For a function periodic on $[-1,1]$, the most natural interpolant is the trigonometric polynomial {eq}`trigcardinalinterp`. However, from {eq}`trigcardinal` one finds that

$$
  :label: trapperiod
  \int_{-1}^1 \sum_{k=-n}^n y_k\tau_k(x)\, dx =  \sum_{k=-n}^n y_k\left[ \int_{-1}^1 \tau_k(x)\, dx\right] = \frac{1}{2n+1} \sum_{k=-n}^n y_k.
$$

In {ref}`prob-specquad-trapperiod` you are asked to verify that this result is identical to the value of the trapezoid formula on $2n+1$ nodes. 
```{index} trapezoid formula
```
 That is, \texthighlight{trapperiod}{the trapezoid integration formula is exponentially accurate for periodic functions.}


::::{prf:example}
  \inputexample{globalapprox}{integrate_ellipse}
::::



## Clenshaw–Curtis integration

Suppose $f$ is smooth but not periodic. If we use a global polynomial interpolating $f$ at the Chebyshev 2nd-kind points from {eq}`chebextreme`,

$$
  %:label: chebextremerepeat
  t_k = - \cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n,
$$

and integrate the resulting polynomial interpolant, the method should have spectral accuracy for smooth integrands. The resulting algorithm is known as 
```{index} Clenshaw–Curtis integration
```
 **Clenshaw–Curtis integration**.

Having specified the nodes in {eq}`globalquad`, all that remains is to find the weights. The Lagrange form of the interpolating polynomial is

$$
p(x) = \sum_{k=0}^{n} f(x_k) \ell_k(x).
$$

From this,
\begin{align*}
  I = \int_{-1}^1 f(x)\, dx \approx \int_{-1}^1 p(x) \d x &= \int_{-1}^1 \sum_{k=0}^n f(x_k) \ell_k(x) \d x\\
    &= \sum_{k=0}^n f(x_k) \int_{-1}^1 \ell_k(x) \d x\\  &= \sum_{k=0}^n c_k f(x_k), \qquad
c_k = \int_{-1}^1 \ell_k(x)\,dx.
\end{align*}
For even values of $n$ the result is

$$
  :label: clencurtweights
  c_k =
\begin{cases}
\dfrac{1}{n^2-1}, & \text{$k=0$ or $k=n$},\\[3mm]
\dfrac{4}{n} \displaystyle \sum_{j=0}^{n/2} \frac{\cos ( 2 \pi j k / n)}{\gamma_j
(1-4 j^2) }, & k=1,\dots,n-1,
\end{cases}
  \quad
\gamma_j =
\begin{cases}
2, & j=0 \text{ or } n/2,\\
1, & j=1,2,\dots,n/2-1.
\end{cases}
$$

There are different formulas for odd values of $n$. Note that all of the weights depend on $n$; e.g.\ $c_2$ for $n=4$ is not the same as $c_2$ for $n=10$. Also note that the interpolant itself never needs to be computed.
%As an example, consider the case with $n=4$.  In
%this case, the weights can be computed with relatively little effort
%with symbolic mathematical software or by hand, and the weights are
%given by
%
$$
%w_k =  \left[ \frac{2}{n+1} +  \frac{2}{n-2+1} \left(
%    \sum_{\substack{j=0\\j\neq k}}^n \
%     \sum_{\substack{m=j+1\\m\ne k}}^n x_j x_m \right) + \frac{2}{n-4+1}
%   \left(  \prod_{\substack{j=0\\j\ne k}}^n x_j \right) \right]
%    \left( \prod_{\substack{j=0\\j\ne k}}^n (x_k - x_j) \right)^{-1},
%$$

%for $k=0,\dots,n$.  This formula is only valid for $n=4$, but
%is written this to suggest the origin of each of the terms: they are all from the
%even terms of $\ell_k$.  The odd terms are integrated over symmetric limits
%and thus drop out.  Substitution for the $x_k$ leads to the rule, for $n=4$,
%
$$
%I \approx \frac{1}{15} \left( y_0 + 8 y_1 + 12 y_2 + 8 y_3 + y_4 \right);
%$$

%however, the value of this approach doesn't lie in such a formula.
{numref}`Function {number}<function-ccint>` performs Clenshaw–Curtis integration for even values of $n$.\footnote{This function is modeled after the function "clencurt.m" of {cite}`TrefSpec`.}

\begin{function}[tbp]
  \caption{("ccint") Clenshaw–Curtis numerical integration.}
  \inputfunction{globalapprox}{ccint}
\end{function}


## Gauss–Legendre integration

Let us reconsider the generic numerical integration formula {eq}`globalquad`,

$$
  \int_{-1}^1 f(x)\, dx \approx \sum_{k=1}^n c_k f(t_k) = Q_{n}[f],
$$

where $Q_n[f]$ stands for the application of the formula to function $f$.
(We now start the sum from $k=1$ instead of $k=0$ for notational convenience in what follows.)
The interpolation approach spurred us to use Chebyshev nodes. But it's far from clear that these are the best nodes for the specific application of finding an integral. Instead, the formula can still be defined as the integral of a polynomial interpolant, but with the nodes and weights all chosen to satisfy an optimality criterion.

The definition of optimality that leads to a very useful result is to require that the formula be exact for all polynomial integrands of as high a degree as possible. Denote the set of all polynomials of degree at most $m$ by $\cp_m$. Since there are $n$ nodes and $n$ weights available to choose, it seems plausible to expect $m=2n-1$, and this turns out to be correct. Hence the goal is now to find nodes $t_k$ and weights $c_k$ such that

$$
  :label: gqoptimality
  \int_{-1}^1 p(x)\,dx = Q_{n}[p] = \sum_{k=1}^n c_k p(t_k), \qquad \text{all $p \in \cp_{2n-1}$}.
$$

If these conditions are satisfied, the resulting method is called 
```{index} Gaussian!integration
```
 **Gauss–Legendre integration**, or often *Gaussian quadrature*.


::::{prf:example}
:label: example-gquad2
  As an example, consider the case $n=2$, which should allow us to satisfy {eq}`gqoptimality` for the polynomials $1$, $x$, $x^2$, and $x^3$. Applying the integration formula in each case, we get the conditions
  
$$
  :label: gausscond2
  \begin{split}
    2 &= c_1 + c_2\\
    0 &= c_1 t_1 + c_2 t_2 \\
    \frac{2}{3} &= c_1 t_1^2 + c_2 t_2^2\\
    0 &= c_1 t_1^3 + c_2 t_2^3.
  \end{split}
  $$

  These equations can be solved to obtain
  
$$
    c_1=c_2=1, \quad x_1 = -\frac{1}{\sqrt{3}}, \quad x_2 =
    \frac{1}{\sqrt{3}},
  $$

  which specifies the two-point Gaussian quadrature formula.
::::


Note that in {ref}`example-gquad2` it was sufficient for the formula to be exact for just the four polynomials $1$, $x$, $x^2$, and $x^3$. Because the formula is linear, i.e., $Q_n[\alpha p + q] = \alpha Q_n[p] + Q_n[q]$, exactness will hold for all of $\cp_3$. This means that the conditions here are sufficient to make $Q_2$ exact for all polynomials in $\cp_3$.

Generalizing the process above to general $n$ would be daunting, as the conditions on the nodes and weights are nonlinear. Fortunately, a more elegant approach is possible.

::::{prf:theorem}
  :label: theorem-gaussquad
  The roots of the Legendre polynomial $P_n(x)$ are the nodes of an $n$-point Gaussian quadrature formula.
::::


::::{prf:proof}
  Choose an arbitrary $p\in\cp_{2n-1}$, and let $I_n[p](x)$ be the interpolating polynomial for $p$ using the as-yet unknown nodes $t_1,\dots,t_n$. By definition,
  \[
    Q_n[p] = \int_{-1}^1 I_n[p](x)\, dx.
  \]
  Since $I_n[p](x)$ has degree $n-1$, it is exactly equal to $p$ if $p\in\cp_{n-1}$, and {eq}`gqoptimality` is trivially satisfied. Otherwise, the error formula {eq}`interperror` implies

$$
  p(x) - I_n[p](x) = \frac{p\spr{n}(\xi(x))}{n!} \Phi(x) = \frac{p\spr{n}(\xi(x))}{n!} (x-t_1)\dots(x-t_n).
$$

Trivially, the left-hand side is a polynomial in $\cp_{2n-1}$ of degree at least $n$, so the right-hand side must be
too. Thus, we can  write

$$
  p(x) - I_n[p](x) = \Psi(x) \Phi(x),
$$

where as above, $\Phi(x)=\prod(x-t_i)$, and $\Psi(x)\in \cp_{n-1}$ is unknown. The optimality requirement {eq}`gqoptimality` becomes

$$
  0 = \int_{-1}^1 p(x)\,dx - Q_{n}[p]  = \int_{-1}^1 \Bigl[p(x) - I_n[p](x)\Bigr]\,dx = \int_{-1}^1 \Psi(x) \Phi(x) \, dx.
$$

Given that $\Psi(x)\in \cp_{n-1}$, we can ensure that this condition is satisfied if

$$
  :label: gqorthogonality
  \int_{-1}^1 q(x)\Phi(x) \,dx = 0 \quad \text{ for all $q \in {\cp}_{n-1}$.}
$$


Hence satisfaction of {eq}`gqorthogonality` implies satisfaction of {eq}`gqoptimality`. But by the orthogonality property of Legendre polynomials, satisfaction of {eq}`gqorthogonality` is guaranteed if $\Phi(x)=cP_n(x)$ for a constant $c$. Thus $\Phi$ and $P_n$ have the same roots.
::::


From {prf:ref}`theorem-legroot` we know that the roots of $P_n$ are distinct and all within $(-1,1)$. (Conversely, it would be strange to have the integral of a function depend on some of its values outside the integration interval!)  There is no explicit formula for the roots. However, they are widely available to high precision, and there are fast algorithms to compute them on demand. {numref}`Function {number}<function-glint>` uses one of the oldest methods for computing the roots and is practical up to a hundred nodes or so.

\begin{function}
  \caption{("glint") Gauss–Legendre numerical integration.}
  \inputfunction{globalapprox}{glint}
\end{function}



## Convergence

Both Clenshaw–Curtis and Gauss–Legendre integration are based on the integration of a global polynomial interpolant, and both are spectrally accurate. The Clenshaw–Curtis method on $n+1$ points is exact for polynomials in $\cp_n$, whereas the Gauss–Legendre method with $n$ points is exact on all of $\cp_{2n-1}$. For this reason, it is possible for Gauss–Legendre to converge at a rate that is "twice as fast," i.e., with roughly the square of the error of Clenshaw–Curtis. But the full story is not simple.


::::{prf:example}
  \inputexample{globalapprox}{quadcompare}
::::


The difference between Clenshaw–Curtis and Gauss–Legendre is often modest—and far less than the difference between spectral and algebraic convergence. It is possible, though, to construct integrands for which adaptivity is critical.  Choosing a method is highly problem-dependent, but a rule of thumb is that for large error tolerances, an adaptive low-order method is likely to be a good choice, while for high accuracy the spectral methods often dominate.

%\fxnote[author=TAD]{Do we want the formulas for Gauss quadrature accuracy?}
% Commented out pending a decision.

%It is possible to analyze the convergence of Clenshaw-Curtis and Gaussian
%quadratures and to see why they are close.  Let $I$ and $I_n$, be defined as, respectively,
%
$$
%I=\int_{-1}f(x)\, dx \text{ and } I_n = \sum_{k=0}^n c_k f(x_k).
%$$

%$I_n$ can be taken to be either the Gaussian or the Clenshaw–Curtis formula.
%For either case, Trefethen (SIAM Review, 50 (2008) 67-87) can compute a bound on the
%error for both methods to be, for $n$ large enough,
%
$$
%|I-I_n| \le \frac{32}{15 \pi} \frac{ ||f^{(j)}||_T }{j(2n+1-j)^j}
%$$

%Here the $T$-norm is the Chebyshev-weighted 1-norm
%\[
%||u||_T = \left| \left| \frac{u'(x)}{\sqrt{1-x^2}} \right| \right|_\infty;
%\]
%for our purposes, it is a measure of the size of the $j$-th derivative in the remainder
%term and it must be finite (bounded).  What is most important here is that the error can
%be analyzed and that the same form for a bound on the error emerges for large $n$.  Trefethen
%explains that this phenomenon can be explained via two different approaches: one is called aliasingand the
%other makes elegant use of the behavior of rational functions and the complex plane.  We refer the
%reader to that article for more details.

% This comparison is moved to the exercises.

% We now follow the approach of Trefethen and compare the behavior of these two methods for a few
%different functions; Figure~\ref{fig:ccgcompconv} shows the results.  All of the functions shown
%are $\int_{-1}^1 f(x)\,dx$.  For $f(x)=x^{24}$, we see that the error plunges to roundoff levels
%at $N=12$ for the Gaussian quadrature, while a similar plunge occurs at $N=24$ for the
%Clenshaw-Curtis method; this illustrates the superior degree of precision for the
%Gaussian approach.  For $f(x)=\cos^2(\pi x)$, there are no singularities even when $x$ is allowed
%to take on complex values ($f$ is said to be **entire** when this is true).  For this $f$, Gaussian
%quadrature converges to the roundoff levels detectably faster than Clenshaw-Curtis.  For $f=1/(1+4x^2)$, which
%we have already seen, the two methods have essentially the same error until
%\begin{figure}
%  \centering
%  \includegraphics[width=6in]{ExaCompCCGConv}
%  \caption{Comparison of Clenshaw–Curtis (solid/circles), Gaussian (dashed/squares) for the integrals of
%  different functions on $[-1,1]$.}
%  \label{fig:ccgcompconv}
%\end{figure}
%about $N=18$, at which point Clenshaw-Curtis takes a slower rate toward
%machine precision.  Trefethen predicts when this turn occurs based on his analysis; the cause is
%related to the fact that there are singlarities in $f$ (at $x=\pm i/2$, $i=\sqrt{-1}$)
%near the interval of integration.
%Finally, $f=|\sin(\pi x)|^3$ has only two continuous derivatives, and the rate of decay in
%each method is equally bad.  This is to be expected with these kinds of methods, which require
%a great deal of smoothness (many derivatives exist everywhere on the integration interval)
%in order to achieve rapid convergence.  For further details, we refer the reader to
%Trefethen's 2008 *SIAM Review* article.


\begin{exercises}
  \input{globalapprox/exercises/SpectralQuadrature}
\end{exercises}

\clearpage
