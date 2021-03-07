# Piecewise linear interpolation

```{index} interpolation; by piecewise polynomials
```

Piecewise linear interpolation is simply a game of connect-the-dots. Let us assume the nodes are given in order, so that $t_0 < t_1 < \cdots < t_n$. Between each pair of adjacent nodes, we use a straight line segment. The resulting interpolant $p(x)$ is given by

```{math}
:label: pwlinear
p(x) = y_k + \frac{y_{k+1}-y_k}{t_{k+1}-t_k}(x-t_k), \quad \text{ for } x\in[t_k,t_{k+1}].
```

It should be clear that on each interval $[t_k,t_{k+1}]$, $p(x)$ is a linear function, and you can easily verify from the formula that it passes through both $(t_k,y_k)$ and $(t_{k+1},y_{k+1})$.

## Hat functions

Rather than basing an implementation on {eq}`pwlinear`, we return to the idea used in {doc}`../linsys/demos/interp-vander` of choosing the interpolant from among the linear combinations of a preselected finite set of functions. In the present context we use

```{math}
:label: hatfun
  H_k(x) =
  \begin{cases}
    \dfrac{x-t_{k-1}}{t_k-t_{k-1}}, & \text{if $x\in[t_{k-1},t_k]$},\\[2.5ex]
    \dfrac{t_{k+1}-x}{t_{k+1}-t_{k}}, & \text{if
      $x\in[t_{k},t_{k+1}]$},\\[2.5ex]
    0, & \text{otherwise},
  \end{cases} \qquad  k=0,\ldots,n.
```

```{index} hat functions
```

The functions $H_0,\ldots,H_n$ are called {term}`hat functions`. They depend on the node vector $\mathbf{t}$, but this dependence is not usually indicated explicitly.

Each hat function is globally continuous and is linear inside every interval $[t_k,t_{k+1}]$.  Consequently, any linear combination of them will have the same property. Furthermore, *any* such function is expressible as a unique linear combination of hat functions, i.e,

```{math}
  :label: plbasis
  \sum_{k=0}^n c_k H_k(x),
```

for some choice of the coefficients $c_0,\ldots,c_n$. No smaller set of functions can have the same properties. We summarize these facts by calling the hat functions a **basis** of the set of functions that are continuous and piecewise linear relative to $\mathbf{t}$.  Another point of view, familiar from abstract linear algebra, is that a basis sets up a one-to-one correspondence between the spanned function space and the more familiar space $\mathbb{R}^{n+1}$, with each function being represented by its coefficients $c_0,\ldots,c_n$.

Note that the definitions of $H_0$ for $x<t_0$ and $H_n$ for $x>t_n$ are irrelevant—a fact exploited by the implementation given in {numref}`Function {number}<function-hatfun>` through the introduction of two fictitious nodes lying on either side of the interval $[t_0,t_n]$. This trick allows the use of an identical formula for all of the hat functions. Otherwise, we would need to take special action for the two edge cases $H_0$ and $H_n$.

(function-hatfun)=

````{proof:function} hatfun
**Hat function/piecewise linear basis function.**

```{code-block} julia
:lineno-start: 1
"""
hatfun(x,t,k)

Evaluate a piecewise linear "hat" function at `x`, where `t` is a
vector of n+1 interpolation nodes and `k` is an integer in 0:n
giving the index of the node where the hat function equals one.
"""

function hatfun(x,t,k)
    n = length(t)-1
    # Return correct node given mathematical index k, including
    # fictitious choices.   
    function node(k)
        if k < 0
            2t[1]-t[2]
        elseif k > n 
            2t[n+1]-t[n] 
        else
            t[k+1]
        end
    end

    H1 = (x-node(k-1))/(node(k)-node(k-1))   # upward slope
    H2 = (node(k+1)-x)/(node(k+1)-node(k))   # downward slope

    H = min(H1,H2)
    return max(0,H)
end
````

````{prf:example} Julia demo
:class: demo
:label: demos-pwlin-hat
{doc}`demos/pwlin-hat`
````

Another trick in {numref}`Function {number}<function-hatfun>` is to exploit the observation that at each $x$, $H_k(x)$ is the larger of two options: zero, or the smaller of two linear functions (those appearing in {eq}`hatfun`). This code is not the most efficient one possible, but it is more concise than detecting which particular subinterval each $x$ lies within.

## Cardinality conditions

```{index} cardinal functions
```

For the purposes of interpolation, the most salient property of the hat functions is that they are cardinal functions for piecewise linear interpolation; that is, they satisfy the **cardinality conditions**

```{math}
:label: cardinalcond
H_k(t_i) =
\begin{cases}
  1, &\text{if $i=k$,}\\
  0, & \text{otherwise.}
\end{cases}
```

The appeal of a cardinal basis is that it makes the expression of the interpolant trivial. All candidate piecewise linear (PL) functions can be expressed as a linear combination such as {eq}`plbasis` for some coefficients $c_0,\ldots,c_n$. But because of the cardinality conditions and the necessity for $p(x)$ to interpolate the data values in $\mathbf{y}$, \texthighlight{cardinalexpress}{expressing the interpolant using the hat functions is trivial:}

```{math}
  :label: plbasissol
  p(x) = \sum_{k=0}^n y_k H_k(x).
```

````{prf:example} Julia demo
:class: demo
:label: demos-pwlin-usage
{doc}`demos/pwlin-usage`
````

The resulting algorithmic simplicity is reflected in {numref}`Function {number}<function-plinterp>`. Take note that the output of {numref}`Function {number}<function-plinterp>` is itself a function, meant to be called with a single argument representing values of $x$. Our mathematical viewpoint is that the result of an interpolation process is a function, and our codes reflect this.

A final appealing characteristic of the hat function basis is that it depends only on the node locations, while the expansion coefficients in {eq}`plbasis` depend only on the data values. This clean separation would be useful if we wanted to construct many interpolants on the same node set, and it has deeper theoretical uses as well.

(function-plinterp)=

````{proof:function} plinterp
**Piecewise linear interpolation.**

```{code-block} julia
:lineno-start: 1
"""
plinterp(t,y)

Create a piecewise linear interpolating function for data values in
`y` given at nodes in `t`.
"""
function plinterp(t,y)
n = length(t)-1
return x -> sum( y[k+1]*hatfun(x,t,k) for k in 0:n )
end
```
````

## Conditioning and convergence

The condition number bounds from \thmref{interp-conditioning} are very simple for piecewise linear interpolation, because the interpolant of the data $\mathbf{e}_k$ is just the hat function $H_k$. Hence $1\le \kappa \le n+1$. However, there is an even simpler result.

```{index} condition number; of interpolation
```

````{prf:theorem} Conditioning of PL interpolation
:label: thm-plcondition
The absolute condition number of piecewise linear interpolation in the infinity norm equals one. More specifically, if $\mathcal{I}$ is the piecewise linear interpolation operator, then 

```{math}
:label: plcondition
\| \mathcal{I}(\mathbf{y}+\mathbf{z}) - \mathcal{I}(\mathbf{y}) \|_\infty = \|\mathbf{z}\|_\infty.
```

(The norm on the left side is on functions, while the norm on the right side is on vectors.)
````

````{prf:proof}
By linearity,

```{math}
\mathcal{I}(\mathbf{y}+\mathbf{z}) - \mathcal{I}(\mathbf{y}) = \mathcal{I}(\mathbf{z}) = \sum_{k=0}^n z_k H_k(x).
```

Call this piecewise linear function $p(x)$. Consider a maximum element of $\mathbf{z}$, i.e. choose $i$ such that $|z_i|=\|\mathbf{z}\|_\infty$. Then $|p(t_i)|=\|\mathbf{z}\|_\infty$. Hence $\|p\|_\infty\ge \|\mathbf{z}\|_\infty$. Now consider

```{math}
|p(x)| = \left|\sum_{k=0}^n z_k H_k(x)\right| \le \sum_{k=0}^n |z_k| H_k(x) \le \|\mathbf{z}\|_\infty \sum_{k=0}^n H_k(x) = \|\mathbf{z}\|_\infty.
```

You are asked to prove the final step above in [an exercise](problem-plpartunity). We conclude that  $\|p\|_\infty\le \|\mathbf{z}\|_\infty$, so that $\|p\|_\infty = \|\mathbf{z}\|_\infty$, which completes the proof.
````

Now suppose that $f$ is a "nice" function on an interval $[a,b]$ containing all of the nodes. We can play a game of sampling values of $f$ to get data, i.e. $y_k=f(t_k)$ for all $k$, then perform piecewise linear interpolation of the data to get a different function, the interpolant $p$. How close is $p$ to the original $f$? To make a simple statement, we will consider only the case of equally spaced nodes covering the interval.

````{prf:theorem} Convergence of PL interpolation
:label: theorem-pl-converge
Suppose that $f(x)$ has a continuous second derivative in $[a,b]$ (often expressed as $f\in C^2[a,b]$). Let $p_n(x)$ be the piecewise linear interpolant of $\bigl(t_i,f(t_i)\bigr)$ for $i=0,\ldots,n$, where $t_i=a+i h$ and $h=(b-a)/n$. Then
  
```{math}
:label: placcuracy
\bigl\| f - p_n \bigr\|_\infty = \max_{x \in [a,b]}
|f(x)-p(x)| \le M h^2,
```

where $M = \bigl\| f'' \bigr\|_\infty$.
````

For an outline of a proof, see [this exercise](problem-placcuracy).

```{margin}
Piecewise linear interpolation is second-order accurate.
```

````{prf:example} Julia demo
:class: demo
:label: demos-pwlin-converge
{doc}`demos/pwlin-converge`
````

We normally don't have access to $f''$, so the importance of {prf:ref}`theorem-pl-converge` is that the error in the interpolant is $O(h^2)$ as $h\to 0$. The leading exponent of 2 is described by saying that piecewise linear interpolation is second-order accurate. For instance, if we double the number of equally spaced nodes used to sample a function, the piecewise linear interpolant becomes about four times more accurate.

## Exercises

1. ⌨ For each of the functions and intervals given, perform piecewise linear interpolation using {numref}`Function {number}<function-plinterp>` for equispaced nodes with $n=10,20,40,80,160$. For each $n$, estimate the error

    ```{math}
    E(n) = \| f-p \|_\infty = \max_x | f(x) - p(x) |
    ```

    by evaluating the function and interpolant at 1600 points in the interval. Make a log--log plot of $E$ as a function of $n$ and add the line $E=Cn^{-2}$ for a constant $C$ of your choosing.

    **(a)** $\cos(\pi^2 x^2)$ on $[0,1]$

    **(b)** $\log(x)$ on $[1,3]$

    **(c)** $\sin(x^2)$ on $[0,2.5]$

    ````{only} solutions
    %%
    %

    %% (a)
    f = @(x) cos((pi*x).^2);
    domain = [0,1];
    xx = linspace(domain(1),domain(2),1600)';
    err = [];
    n_ = 10*[1 2 4 8 16];
    for n = n_
        h = diff(domain)/n;
        x = domain(1) + h*(0:n)';
        p = plinterp(x,f(x));
        err = [ err; norm(p(xx)-f(xx),inf) ];
    end
    fprintf('   n      err\n')
    fprintf(' %3i   %6.2e\n',[ n_; err' ]);
    fprintf('\n\n')
    clf, loglog(n_,err,'-o')
    hold on, loglog(n_,60*n_.^(-2),':')
    axis tight, xlabel('n'), ylabel('e_n')
    title(['function ' char(f)])

    %% (b)
    f = @log;
    domain = [1,3];
    xx = linspace(domain(1),domain(2),1600)';
    err = [];
    n_ = 10*[1 2 4 8 16];
    for n = n_
        h = diff(domain)/n;
        x = domain(1) + h*(0:n)';
        yy = plinterp(xx,x,f(x));
        err = [ err; norm(yy-f(xx),inf) ];
    end
    fprintf('   n      err\n')
    fprintf(' %3i   %6.2e\n',[ n_; err' ]);
    fprintf('\n\n')
    clf, loglog(n_,err,'-o')
    hold on, loglog(n_,60*n_.^(-2),':')
    axis tight, xlabel('n'), ylabel('e_n')

    %% (c)
    f = @(x) sin(x.^2);
    domain = [0,2.5];
    xx = linspace(domain(1),domain(2),1600)';
    err = [];
    n_ = 10*[1 2 4 8 16];
    for n = n_
        h = diff(domain)/n;
        x = domain(1) + h*(0:n)';
        yy = plinterp(xx,x,f(x));
        err = [ err; norm(yy-f(xx),inf) ];
    end
    fprintf('   n      err\n')
    fprintf(' %3i   %6.2e\n',[ n_; err' ]);
    fprintf('\n\n')
    clf, loglog(n_,err,'-o')
    hold on, loglog(n_,60*n_.^(-2),':')
    axis tight, xlabel('n'), ylabel('e_n')
    ````

2. ✍ For this problem, let $H(x)$ be the hat  function that passes through the three points $(-1,0)$, $(0,1)$, and $(1,0)$.

    **(a)** Write out a piecewise definition of $H$ in the style of {eq}`hatfun`.

    **(b)** Define the function $Q$ by $Q(x) = \int_{x-1}^x H(t)\, dt$. Find a piecewise formula for $Q(x)$. (Hint: Perform the integration separately for the cases $-1\le x \le 0$, $0\le x \le 1$, etc.)

    **(c)** Make a sketch of $Q(x)$ for $-2\le x \le 2$.

    **(d)** Show that $Q$ is continuous. Are $Q'$ and $Q''$?
  
    ````{only} solutions
    % Over the three subintervals $[-1,0]$, $[0,1]$, $[1,2]$, the results are $\tfrac{1}{2}+x+\tfrac{1}{2}x^2$, $\tfrac{1}{2}+x-x^2$, and $2-2x+\tfrac{1}{2}x^2$. $Q$ and $Q'$ are continuous, but $Q''$ is not.
    ````

3. ✍ Before electronic calculators, the function $\ln(x)$ was often computed using piecewise linear interpolation with a table of values. If you were using such a table at the nodes $1.01,1.02,\ldots,1.99,2$, what is an upper bound on the error in the result?

    ````{only} solutions
    % Use Theorem 5.2.1 to bound the error with $M=1$. Since $h=10^{-2}$, the error is bounded by $10^{-4}$.
    ````

    (problem-plpartunity)=
4. ✍ Show that for any node distribution and any $x\in[t_0,t_n]$,
  
    ```{math}
    :label: plpu
    \sum_{k=0}^n H_k(x) = 1.
    ```

    (Hint: The simplest way is to apply {eq}`plbasissol`.) This is called the **partition of unity** property.

    ````{only} solutions
    ````

    (problem-placcuracy)=
5. ✍ Here we consider a proof of the {prf:ref}`theorem-pl-converge` using the mean value theorems from elementary calculus: If $f$ is continuously differentiable in $(a,b)$, then there exist points $s$ and $t$ in $(a,b)$ such that
  
    ```{math}
    \int_a^b f(z) \, dz = (b-a)f(s) \qquad \text{and} \qquad f'(t) = \frac{f(b)-f(a)}{b-a}.
    ```

    For the following, suppose $x \in (t_k,t_{k+1})$.

    **(a)** Show that for some $s \in (t_k,t_{k+1})$,

    ```{math}
    f(x) = y_k + (x-t_k)f'(s).
    ```

    **(b)** Show that for some other values $u$ and $v$ in $(t_k,t_{k+1})$,

    ```{math}
    f'(s) -  \frac{y_{k+1}-y_k}{t_{k+1}-t_k} = (s-u) f''(v).
    ```

    **(c)** Use {eq}`pwlinear` to finish the proof of the theorem.
  
