# The barycentric formula

The Lagrange formula {eq}`lagrangeinterp` is useful theoretically but not ideal for computation. For each new value of $x$, all of the cardinal functions $\ell_k$ must be evaluated at $x$, which requires a product of $n$ terms. Thus the total work is $O(n^2)$ for every value of $x$. Moreover, the formula is numerically unstable (see [this exercise](problem-lagrange-instability)). An alternative version of the formula improves on both issues.

## Derivation

Define, as in the proof of {prf:ref}`theorem-interperror`,

:::{math}
  :label: productpoly
  \Phi(x) = \prod_{j=0}^n (x-t_j).
:::

```{index} barycentric weights
```
Also define the **barycentric weights**

:::{math}
  :label: baryweight
  w_k = \frac{1}{\displaystyle \prod_{\substack{j=0\\j\neq k}}^n (t_k - t_j)}, \qquad
  k = 0,\ldots,n.
:::

Then {eq}`lagrange` can be written as

:::{math}
  :label: lagrangealt
  \ell_k(x) = \Phi(x) \frac{w_k}{x-t_k},
:::

and thus the interpolating polynomial for $f(x)$ is

:::{math}
  :label: bary1
  p(x) = \Phi(x) \sum_{k=0}^n \frac{w_k}{x-t_k} y_k.
:::

We are still one step away from the most useful formula. Obviously, the constant function $p(x)\equiv 1$ is its own polynomial interpolant on any set of nodes. The uniqueness of the interpolating polynomial, as proved in {prf:ref}`theorem-polyinterp`, allows us to plug $y_k=1$ for all $k$ into {eq}`bary1` to obtain

$$
  1 = \Phi(x) \sum_{k=0}^n \frac{w_k}{x-t_k}.
$$

This is solved for $\Phi(x)$ and put back into {eq}`bary1` to get

:::{math}
  :label: bary2
  p(x) = \frac{\displaystyle \sum_{k=0}^n y_k \, \dfrac{w_k}{x-t_k}  }{\displaystyle\sum_{k=0}^n \dfrac{w_k}{x-t_k}},
:::

```{index} barycentric interpolation formula
```
which is the {term}`barycentric formula` for the interpolating polynomial.[^bary12] Equation {eq}`bary2` is certainly an odd-looking way to write a polynomial! But the barycentric formula is the key to efficient and stable evaluation of a polynomial interpolant.

[^bary12]: More precisely,  {eq}`bary1` and  {eq}`bary2` are called the *first* and *second* barycentric formulas, respectively.

::::{prf:example}
:label: example-writeoutbary2
Let us write out the barycentric formula for the interpolating polynomial for the quadratic case ($n=2$) for {prf:ref}`example-ClassicalLagrange`.  The weights are computed from {eq}`baryweight`:
  
:::{math}
  w_0 = \frac{1}{(t_0-t_1)(t_0-t_2)} = \frac{1}{\left(0-\frac{\pi}{6}\right)
\left(0-\frac{\pi}{3}\right)} = \frac{18}{\pi^2},
:::

and similarly, $w_1 = -36/\pi^2$ and $w_2=18/\pi^2$.

Note that in {eq}`bary2`, any common factor in the weights cancels out without affecting the results. Hence it's a lot easier to use $w_0=w_2=1$ and $w_1=-2$. Then

\begin{align*}
    p(x) & = \frac{\rule[-1.2em]{0pt}{1em} \dfrac{w_0}{x-t_0} y_0  + \dfrac{w_1}{x-t_1} y_1 + \dfrac{w_2}{x-t_2} y_2 }{ \rule{0pt}{1.5em} \dfrac{w_0}{x-t_0} + \dfrac{w_1}{x-t_1} + \dfrac{w_2}{x-t_2}}\\[1.5ex]
    & =\frac{ \rule[-1.2em]{0pt}{1em}\left( \dfrac{1}{x} \right) 0 -  \left( \dfrac{2}{x-\pi/6} \right) \dfrac{1}{\sqrt{3}} + \left( \dfrac{1}{x-\pi/3} \right) \sqrt{3} }{
        \rule{0pt}{1.6em} \dfrac{1}{x} - \dfrac{2}{x-\pi/6} + \dfrac{1}{x-\pi/3}  }
\end{align*}
  
Further algebraic manipulation could return this expression to the classical Lagrange form derived in {prf:ref}`example-ClassicalLagrange`.
::::

For certain canonical node distributions, simple formulas for the weights $w_k$ are known. Otherwise, computing all $n+1$ weights from {eq}`baryweight` takes $O(n^2)$ operations. However, the weights depend only on the nodes, not the data—and once they are known, computing $p(x)$ from {eq}`bary1` for any set of data at a particular value of $x$ takes just $O(n)$ operations.

## Implementation

(function-polyinterp)=
````{proof:function} polyinterp
**Polynomial interpolation by the barycentric formula**

```{code-block} julia
:lineno-start: 1
"""
polyinterp(t,y)

Return a callable polynomial interpolant through the points in
vectors `t`,`y`. Uses the barycentric interpolation formula.
"""
function polyinterp(t,y)
    @assert (isa(t,OffsetArray) && isa(y,OffsetArray)) "Vectors must be indexed 0:n"
    n = length(t)-1
    C = (t[n]-t[0]) / 4           # scaling factor to ensure stability
    tc = t/C
    
    # Adding one node at a time, compute inverses of the weights.
    ω = OffsetArray(ones(n+1),0:n)
    for m in 0:n-1
        d = tc[0:m] .- tc[m+1]    # vector of node differences
        @. ω[0:m] *= d            # update previous
        ω[m+1] = prod( -d )       # compute the new one
    end
    w = 1 ./ ω                    # go from inverses to weights

    p = function (x)
        # Compute interpolant.
        terms = @. w / (x - t)
        if any(isinf.(terms))     # there was division by zero
            # Apply L'Hôpital's Rule exactly.
            idx = findfirst(x.==t)
            f = y[idx]
        else
            f = sum(y.*terms) / sum(terms)
        end
    end
    return p
end
```
````

In {numref}`Function {number}<function-polyinterp>` we show an implementation of the barycentric formula for polynomial interpolation. The first phase is to compute the weights $w_k$, or more conveniently, $w_k^{-1}$. As noted in {prf:ref}`example-writeoutbary2`, a common scaling factor in the weights does not affect the barycentric formula {eq}`bary2`. In our code this fact is used to rescale the nodes so as to avoid arriving at very small or very large numbers.

The weight computation begins with the singleton node set $\{t_0\}$, for which one gets the single weight $w_0=1$. The idea is to grow this single node into the set of all the nodes through a recursive formula. Define $\omega_{k,m-1}$ (for $k< m$) as the inverse of the weight for node $k$ using the set $\{t_0,\ldots,t_{m-1}\}$. Then

$$
  \omega_{k,m} = \displaystyle \prod_{\substack{j=0\\j\neq k}}^{m} (t_k - t_j)
     = \omega_{k,m-1} \cdot (t_k-t_{m}), \qquad k=0,1,\ldots,m-1.
$$

A direct application of {eq}`baryweight` can be used to find $\omega_{m,m}$. This process is iterated over $m=1,\ldots,n$ to find $w_k=\omega_{k,n}^{-1}$.

Once the weights are computed, the function loops over the interpolation nodes to compute the sums in the numerator and denominator of formula {eq}`bary2`. Finally, the code addresses a peculiar feature of {eq}`bary2`: if $x=t_i$ for some value of $i$, the formula calls for division by zero.  Analytically, L'Hôpital's rule applies, and the interpolant takes the prescribed data value for node $t_i$ (see [this exercise](problem-comppoly-barylimit)). Because division by zero numerically causes a `NaN` (not a number) value to be produced, however, we look for such cases and revise them accordingly.

::::{prf:example} Julia demo
:class: demo
:label: demos-barycentric-example
{doc}`demos/barycentric-example`
::::

## Stability

You might expect that as the evaluation point $x$ approaches a node $t_i$, subtractive cancellation error will creep into the barycentric formula because of the term $1/(x-t_i)$. While such errors do occur, they turn out not to cause trouble, though, because the *same* cancellation happens in the numerator and denominator. In fact the barycentric formula is a stable way to evaluate the interpolating polynomial, a statement which we do not try to prove here.

<!-- 
\begin{exercises}
	\input{globalfuncapprox/exercises/Barycentric}
\end{exercises}
 -->
