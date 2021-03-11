# Stability of polynomial interpolation
\label{sec:runge}


```{index} stability!of polynomial interpolation
```

With  barycentric interpolation available in the form of {numref}`Function {number}<function-polyinterp>`, we can explore polynomial interpolation using a numerically stable algorithm. Any remaining sensitivity to error is due to the interpolation process itself.

::::{prf:example}
  \inputexample{globalapprox}{equiinterp1}
::::


%Similar results are obtained for the function from {ref}`example-vander2`. This time, we also investigate the error as a function of $x$.
%
%
::::{prf:example}
%  \inputexample{globalapprox}{equiinterp2}
%::::


%It's useful to compare the results of {ref}`example-vander2` and {ref}`example-equiinterp2`, which are based on interpolation of the same function. When we increased $n$ in the Vandermonde approach, the errors decreased to about $10^{-2}$ and then slowly grew. Using the barycentric formula, the errors did not grow.
%The figure in that example suggests that polynomial interpolants do have the capability of converging, at least in the central part of the domain. So instead we will investigate the possibility of changing the locations of the nodes. Since the regions near the boundaries seem to cause the trouble, it stands to reason that we want to concentrate more nodes there.



## Runge phenomenon

The disappointing loss of convergence in {ref}`example-equiinterp1` is due to the use of equally spaced nodes. We will examine this effect using the error formula {eq}`interperror` as a guide:

$$
  f(x) - p(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!} \Phi(x), \qquad \Phi(x) =
  \prod_{i=0}^n (x-t_i).
$$

The $\Phi(x)$ term can be studied as a function of the nodes only.


::::{prf:example}
  \inputexample{globalapprox}{equiPhi}
::::


Even though $\Phi\to 0$ at every point in the interval, the exponentially growing gap between the ends and the middle of the interval can ruin the convergence of polynomial interpolation for many choices of $f$.


::::{prf:example}
  \inputexample{globalapprox}{Runge}
::::


The observation of instability in {ref}`example-Runge` is known as the 
```{index} Runge phenomenon
```
 {term}`Runge phenomenon`. \texthighlight{runge}{The Runge phenomenon is an instability in the abstract mapping from a function to its polynomial interpolant, manifested when the nodes of the interpolant are equally spaced and the degree of the polynomial increases.}  We reiterate that the phenomenon is rooted in the convergence theory and not a consequence of the algorithm chosen to implement polynomial interpolation.

Significantly, the convergence observed in {ref}`example-Runge` is stable within a middle portion of the interval. By redistributing the interpolation nodes, we will next sacrifice a little of the convergence in the middle portion in order to improve it enough near the ends to rescue the process globally.

## Chebyshev nodes

The observations above suggest that we might find success by having more nodes near the ends of the interval than in the middle. Though we will not give the details, it turns out that there is a precise asymptotic sense in which this must be done to make polynomial interpolation work over the entire interval. \texthighlight{chebpoints}{One especially important node family that gives stable convergence for polynomial interpolation is the 
```{index} Chebyshev points!second kind
```
 **Chebyshev points of the second kind**} (also known as Chebyshev extreme points) defined by

$$
  :label: chebextreme
  t_k = - \cos\left(\frac{k \pi}{n}\right), \qquad k=0,\ldots,n.
$$

These are the projections onto the $x$-axis of $n$ equally spaced points on a unit circle. As such they are densely clustered near the ends of $[-1,1]$, and this turns
out to overcome the Runge phenomenon.

%In MATLAB one can create a vector of these nodes in just one line:
%\begin{verbatim}
%x = -cos( (0:n)'*pi/n );
%\end{verbatim}
%If we repeat the experiment of {ref}`example-equiinterp` with Chebyshev points as the interpolation nodes, the results are very different.


::::{prf:example}
  We repeat {ref}`example-equiPhi` but replace equally spaced nodes with Chebyshev points.
  \inputexample{globalapprox}{chebPhi}
::::



::::{prf:example}
  We repeat {ref}`example-Runge`, replacing equally spaced nodes with Chebyshev nodes.
  \inputexample{globalapprox}{chebRunge}
::::


As a bonus, for Chebyshev nodes the barycentric weights are simple:

$$
  :label: weightcheb
  w_k = (-1)^k d_k, \qquad d_k =
  \begin{cases}
    1/2, & \text{if $k=0$ or $k=n$},\\
    1, & \text{otherwise}.
  \end{cases}
$$




## Spectral convergence
\texthighlight{spectralinterp}{If we take $n\rightarrow \infty$ and use polynomial interpolation on Chebyshev nodes, the convergence rate is exponential in $n$.} The following \fxnote{Credit SMM.} is typical of the results that can be proved.

::::{prf:theorem}
  :label: theorem-spectral
  Suppose $f(x)$ is analytic in an open real interval containing $[-1,1]$. Then there exist constants $C>0$ and $K>1$ such that
  
$$
    :label: spectral
    \max_{x\in[-1,1]} | f(x) - p(x) | \le C K^{-n},
  $$

  where $p$ is the unique polynomial of degree $n$ or less defined by interpolation on $n+1$ Chebyshev 2nd kind points.
::::

The condition "$f$ is analytic" means that the Taylor series of $f$ converges to $f(x)$ in an open interval containing $[-1,1]$.\footnote{Alternatively it means the function is extensible to one that is differentiable in the complex plane.} A necessary condition of analyticity is that $f$ is infinitely differentiable.

In some contexts we refer to {eq}`spectral` as linear convergence, but here it is typical to say that the rate is exponential, geometric, or 
```{index} convergence rate!spectral
```
 {term}`spectral convergence`. One achieves  constant reduction factors in the error by constant increments of $n$. By contrast, algebraic convergence in the form $O(n^{-p})$ for some $p>0$ requires *multiplying* $n$ by a constant factor in order to reduce error by a constant factor. Graphically, spectral error is linear on a log–linear scale, while algebraic convergence is a straight line on a log–log scale.

\begin{exercises}
  \input{globalapprox/exercises/NodeLocations}
\end{exercises}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\clearpage
% # Chebyshev polynomials
% \label{sec:cheb}
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** Three related variables
% *** Definition and formulas
% *** (foreshadowing) FFT
% *** Convergence properties

% ## Real interval and unit circle

% Some of the most fascinating facts about Chebyshev polynomials are
% exposed through the change of variable $x=\cos \theta$. For $-1\le x
% \le 1$, which is where the original polynomial interpolation takes
% place, we can think of $0\le \theta \le 2\pi$ as the angle around the
% unit circle, which is periodic. Both the first-kind and second-kind Chebyshev nodes,
% which seem to have complicated spacing on the interval $[-1,1]$, are
% just equally spaced in $\theta$, as illustrated in
% \begin{figure}
%   \centering
%   \includegraphics{chebvars}
%   \caption{Unit interval and unit circle for Chebyshev polynomials. The interval $[-1,1]$ is the domain of the original interpolation variable $x$. The second-kind points (dots) and the first-kind points (stars) are equally spaced in $\theta$, the angle around the circle, but crowded near the endpoints in $x$.}
%   \label{fig:chebvars}
% \end{figure}

% Referring back to equation~(\ref{eq:chebpoly}), it's clear that when $T_n(x)$ is expressed in $\theta$, it's just $\cos(n\theta)$. Because we can express the degree-$n$ polynomial interpolant $p(x)$ using the Chebyshev basis, we can write
% 
$$
%   p(x) = \sum_{k=0}^n c_k T_k(x) = \sum_{k=0}^n c_k \cos (k \theta).
% $$

% This last sum is a particular type of *trigonometric polynomial*. As we will see in section~\ref{sec:triginterp} (or as you may know from Fourier series), we usually see terms of the form $\sin(k\theta)$ in these as well. However, because each point on the interval $[-1,1]$ maps to two symmetrically located points on the circle, symmetry implies that the cosine terms are sufficient.



% ## Exercises
% \begin{exercises}


% \end{exercises}


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% \clearpage

\clearpage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
