#  Improper integrals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\label{sec:improper}

% \fxnote{Singularity section still needs some work.}

When the interval of integration or the integrand function is unbounded, we say an integral is *improper*. Improper integrals present particular challenges to numerical computation.

In the case of an infinite interval, infinitely many nodes are needed to represent the integrand everywhere, which is absurd. However, in order for the integral to be finite, the integrand has to decay at infinity, which brings up the possibility of truncating the interval:

$$
  :label: infiniteinterval
  \int_{-\infty}^{\infty} f(x)\, dx \approx  \int_{-M}^{M} f(x)\, dx.
$$

This integral can be discretized finitely by, say, the trapezoid formula. Yet this approach might not be realistic.


::::{prf:example}
:label: example-slowdecay
  Consider the case with $f(x)=1/(1+x^2)$,
  
$$
    \int_{-\infty}^\infty \frac{1}{1+x^{2}}\, dx = \pi.
  $$

  Although $f$ decays fast enough for the integral to be finite, the decay is not fast enough for simple truncation. Note that $\int_M^\infty f\,dx\approx \int_M^\infty x^{-2}\,dx =  M^{-1}$ when $M$ is large compared to 1. To get 6 digits of accuracy, then, we need to truncate with $M>10^6$. A trapezoid discretization of {eq}`infiniteinterval` with $h<1$ would need millions of nodes.
::::


We can improve the situation a great deal by introducing a variable transformation $x(t)$. In the new variable $t$, the integrand becomes $f(x(t))x'(t)$, and the combination can be made to decay much more rapidly than $f(x)$ does.


## Review of hyperbolic functions

The variable transformation we will use is most naturally stated in terms of the hyperbolic functions, so it's worth a quick refresher on those. Recall the definitions

$$
  :label: hyperdef
  \sinh(t) = \frac{e^t-e^{-t}}{2}, \quad
  \cosh(t) = \frac{e^t+e^{-t}}{2}.
$$

Away from $t\approx 0$, these functions behave essentially like exponentials:

$$
  :label: hyperasymp
  \sinh(t) \approx \pm \frac{1}{2} e^{|t|}, \quad
  \cosh(t) \approx \frac{1}{2} e^{|t|}, \quad\:	
  \text{as } t\to\pm\infty.
$$

The identities

$$
  :label: hyperident
  \frac{d }{d t} \sinh(t) = \cosh(t), \quad \frac{d }{d t} \cosh(t) = \sinh(t), \quad \cosh^2(t)-\sinh^2(t)=1
$$

frequently are handy. Finally, we recall the hyperbolic tangent,

$$
  :label: hypertan
  \tanh(t) = \frac{\sinh(t)}{\cosh(t)}\to\pm 1, \qquad \text{as } t\to\pm\infty.
$$


## Doubly exponential transformation

Returning to the integral {eq}`infiniteinterval`, a particularly useful way to change the integration variable is

$$
  :label: DEquadtrans1
  x = \sinh\left( \frac{\pi}{2} \sinh t \right).
$$

Note that $x=0$ when $t=0$, and $x\rightarrow\pm\infty$ as
$t\rightarrow\pm\infty$. More specifically,

$$
  x \approx \pm \frac{1}{2}e^{\frac{\pi}{4} e^{|t|}}, \qquad \text{as } t\rightarrow\pm\infty,
$$

and {eq}`DEquadtrans1` is often referred to as a 
```{index} doubly exponential transformation
```
 **doubly exponential** transformation.

By the chain rule,

$$
  :label: DEquadchain1
  \int_{-\infty}^\infty f(x)\, dx = \int_{-\infty}^\infty f(x(t))\frac{dx}{dt}\, dt
  = \frac{\pi}{2}\int_{-\infty}^\infty f(x(t))\,
  \cosh\left( \frac{\pi}{2} \sinh t \right)  \cosh t  \,
  dt.
$$

The exponential terms introduced by the chain rule grow doubly exponentially, so we seem to have made the integral much more difficult! But the decay of $f$ in the new variable more than makes up for the new terms, and \texthighlight{detrunc}{doubly exponential transformation makes truncating an infinite interval easy.}


::::{prf:example}
  Consider again the case $f(x)=1/(1+x^2)$ from {ref}`example-slowdecay`. Suppose $x=x(t)$ as in {eq}`DEquadtrans1`. Although the chain rule term is doubly exponential in $t$, $x$ itself is doubly exponential, and since it's squared and in the denominator, the integrand in {eq}`DEquadchain1` *decays* doubly exponentially.
  \inputexample{globalfuncapprox}{dedecay}
::::


We apply the trapezoid rule to the truncated integral in the new variable $t$:

$$
  :label: unboundtrap
  \int_{-\infty}^{\infty} f(x)\, dx = \int_{-\infty}^{\infty} g(t)\, dt \approx h \sum_{k=-K}^K g(kh),
$$

where $g$ consists of $f$ and the chain rule terms from {eq}`DEquadchain1`. This implies that we truncate the integral in $t$ at $t=\pm Kh$. If we interpret this in terms of the original variable, the truncation point is $\pm M=\pm x(Kh)$. As $M\to\infty$ we get the rule of thumb

$$
  :label: DEfindk
  M \approx \frac{1}{2}e^{\frac{\pi}{4} e^{Kh}}, \text{ or }
  K \approx \frac{1}{h}\log\left(\frac{4}{\pi} \log 2M \right).
$$

{numref}`Function {number}<function-intde>` implements doubly exponential integration using $h$ and $M$ as input parameters of the discretization.  The 
{term}`ceil` function helps determine the truncation point.

\begin{function}[tbp]
  \caption{("intde") Doubly exponential integration over $(-\infty,\infty)$.}
  \inputfunction{globalfuncapprox}{intde}
\end{function}



::::{prf:example}
  We use doubly exponential transformation to integrate $\pi=\int_{-\infty}^\infty 1/(1+x^2)\, dx$. We let $h$ and $M$ vary over orders of magnitude and measure the effect on the accuracy of the result from {numref}`Function {number}<function-intde>`.
  \inputexample{globalfuncapprox}{unboundint}
::::



## Integrand singularities

If $f$ asymptotically approaches infinity as $x$ approaches an integration interval endpoint, its exact integral may or may not be finite. Even if $f$ is integrable, however, the methods we have used so far will generally not work well.
%The convergence situation may not be much better if $f$ is finite but has a derivative that is unbounded at an endpoint.


::::{prf:example}
  Let's use {numref}`Function {number}<function-intadapt>` to try to integrate the function $\sqrt{x}$ over $[0,1]$. The exact result is 2/3.
  \inputexample{globalfuncapprox}{intsqrtadapt}
::::


Doubly exponential transformations can be helpful for integrals with such singularities. The change of variable

$$
  :label: DEquadtrans2
  x = \tanh\left( \frac{\pi}{2} \sinh t \right)
$$

is a transformation between $x\in(-1,1)$ and $t\in(-\infty,\infty)$. One finds now that

$$
  :label: DEquadchain2
  \int_{-1}^1 f(x)\, dx = \int_{-\infty}^\infty f(t)\frac{dx}{dt}\, dt
  = \frac{\pi}{2}\int_{-\infty}^\infty f(t)
  \cosh t  \left[\cosh\left( \frac{\pi}{2} \sinh t \right)\right]^{-2}\, dt.
$$

Now there is a squared doubly exponential term in the denominator, so that \texthighlight{desing}{the doubly exponential transformation is able to cancel out even unbounded growth in $f$ as one approaches an endpoint.}

Again in the new $t$ variable we have an unbounded domain that is to be truncated before applying the trapezoid rule, as in {eq}`unboundtrap`, defined by the discretization size $h$ and the upper sum limit $K$.
Since the parameter $K$ is not simple to choose, we will swap it for another one that is easier to understand. Let us truncate
the original integration interval in $x$ to $[-1+\delta,1-\delta]$ for $\delta\ll 1$. One can show that for large positive values of $s$,

$$
  :label: tanhlarge
  \tanh(s) \approx 1 - 2 e^{-2s},
$$

and hence the truncated endpoints in $x$ map to
$\log\bigl[-(2/\pi)\log(\delta/2)\bigr]$ in $t$. If this equals
the rightmost trapezoid node at $t=hK$, we have

$$
  K \approx \frac{1}{h} \log \left[ -\frac{2}{\pi} \log\left(\frac{\delta}{2}\right) \right].
$$

This formula is used in {numref}`Function {number}<function-intsing>`, which accepts $h$ and $\delta$ to define the discretization.
\begin{function}[tbp]
  \caption{("intsing") Integrate function with endpoint singularities.}
  \inputfunction{globalfuncapprox}{intsing}
\end{function}


::::{prf:example}
  \inputexample{globalfuncapprox}{intsqrtde}
::::


Doubly-exponential mapped integration is a remarkably effective technique for a variety of integrals on unbounded domains or with singular integrands. It can also be beneficial for easier problems. Choosing the parameters well can be a nontrivial matter in the general case, however.


\begin{exercises}
  \input{globalfuncapprox/exercises/DoublyExp}
\end{exercises}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\clearpage
\subsection*{Key ideas in this chapter}
\begin{remunerate}
\item There is a unique constructible polynomial of degree less than $n+1$ interpolating $n+1$ points (\hiref{polyexistunique})
\item The Lagrange cardinal polynomials give a simple (though unstable) expression for the interpolating polynomial (\hiref{lagrange}).
\item There is a useful formula for the error in a polynomial interpolant, when the data are samples of a smooth function (\hiref{interperror}).
\item The barycentric formula is the key to efficient and stable evaluation of a polynomial interpolant (\hiref{barycentric}).
\item The Runge phenomenon is an instability in the abstract mapping from a function to its polynomial interpolant, manifested when the nodes of the interpolant are equally spaced and the degree of the polynomial increases (\hiref{runge}).
\item One especially important node family that gives stable convergence for polynomial interpolation is the Chebyshev points of the 2nd kind (\hiref{chebpoints}).
\item If we let the degree $n\rightarrow \infty$ and use polynomial interpolation on Chebyshev nodes, the convergence rate is exponential in $n$ (\hiref{spectralinterp}).
\item We can extend least-squares fitting from data to functions by extending several familiar finite-dimensional definitions to functions (\hiref{continuousextend}).
\item The Legendre polynomials are orthogonal with respect to function inner products (\hiref{legendreorth}).
\item The Chebyshev polynomials are orthogonal with respect to an inner product having a non-unit weight function (\hiref{cheborth}).
\item Trigonometric interpolation can use equally spaced nodes without any difficulty (\hiref{trigequi}).
\item The FFT solves for $N$ trigonometric interpolation coefficients in $O(N\log N)$ flops (\hiref{fft}).
\item When an integrand is interpolated by a spectrally accurate method, the resulting integration formula (e.g., Clenshaw–Curtis or Gauss–Legendre) is also spectrally accurate (\hiref{specint}).
\item The trapezoid integration formula is exponentially accurate for periodic functions (\hiref{trapperiod}).
\item A doubly exponential variable transformation makes truncating an infinite integration interval easy (\hiref{detrunc}).
\item A doubly exponential transformation is able to cancel out a singularity in $f$ as one approaches an integral endpoint (\hiref{desing}).
\end{remunerate}

\subsection*{Where to learn more}

