# Efficiency of matrix computations

Predicting how long an algorithm will take to solve a particular problem, on a particular computer, as written in a particular way in a particular programming language, is an enormously difficult undertaking. It's more practical to predict how the required time will scale as a function of the size of the problem. In the case of a linear system of equations, the problem size is $n$, the number of equations variables in the system.  Because expressions of computational time are necessarily approximate, it's customary to suppress all but the term that is dominant as $n\to\infty$. We first need to build some terminology for these expressions.

```{index} asymptotic notation
```

## Asymptotic notation

For positive-valued functions $f(n)$ and $g(n)$, we say $f(n)=O\bigl(g(n)\bigr)$ ("$f$ is {term}`big-O` of $g$") as $n\rightarrow \infty$ if $f(n)/g(n)$ is bounded above for all sufficiently large $n$. We say $f(n)\sim g(n)$ ("$f$ is {term}`asymptotic` to $g$") as $n\rightarrow \infty$ if $f(x)/g(x)\rightarrow 1$ as $n\rightarrow\infty$. One immediate result is that $f\sim g$ implies $f=O(g)$.[^sets]

[^sets]: More precisely, $O(g)$ and $\sim g$ are *sets* of functions, and $\sim g$ is a subset of $O(g)$. That we write $f=O(g)$ rather than $f\in O(g)$ is a quirk of convention.

````{prf:example}
Consider the functions $f(n) = a_1 n^3 + b_1 n^2 + c_1 n$ and $g(n) = a_2 n^3$ in the limit $n\to \infty$.  Then
  
```{math}
    \lim_{n \to \infty} \frac{f(n)}{g(n)}
    = \lim_{n \to \infty} \frac{a_1 + b_1n^{-1} + c_1n^{-2}}{a_2} =
    \frac{a_1}{a_2} .
```

Since $a_1/a_2$ is a constant, $f(n) = O(g(n))$; if $a_1=a_2$, then $f \sim g$.
````

````{prf:example}
  Consider $f(n) = \sin (1/n)$, $g(n)=1/n$ and $h(n) = 1/n^2$. For large $n$, Taylor's theorem with remainder implies that
  
```{math}
f(n) = \frac{1}{n} - \cos(1/\xi)\frac{1}{6 n^3},
```

where $n<\xi<\infty$.  But
  
```{math}
\lim_{n\to \infty} \frac{f}{g} = \lim_{n\to \infty} 1-\cos(1/\xi)\frac{1}{6 n^2} = 1,
```

and so $f \sim g$.  On the other hand, comparing $f$ and $h$, we find

```{math}
\lim_{n\to \infty} \frac{f}{h} = \lim_{n\to \infty}  n-\cos(1/\xi)\frac{1}{6 n} = \infty,
```

so we cannot say that $f = O(h)$. A consideration of $h/f$ will show that $h = O(f)$, however.
````

It's conventional to use asymptotic notation that is as specific as possible. For instance, while it is true that $n^2+n=O(n^{10})$, it's more informative, and usually expected, to say $n^2+n=O(n^2)$. There are additional notations that enforce this requirement strictly, but we will just stick to the informal understanding.

## Flop counting

```{index} flops
```

Traditionally, in numerical linear algebra we count **floating point operations**, or {term}`flops` for short. In our interpretation each scalar addition, subtraction, multiplication, division, and square root counts as one flop. Given any algorithm, we can simply add up the number of scalar flops and ignore everything else.

```{prf:example} Julia demo
:class: demo
{doc}`demos/flops-mvmult`

{doc}`demos/flops-loglog`
```

It's clear from the example in {doc}`demos/flops-mvmult` that the runtime increases at a function of $n$—but at what rate? Suppose that the time obeys a function that is not just $O(n^p)$, but actually equal to $Cn^p$ for some constants $C$ and $p$. For large enough $n$, this should be a good approximation. Then

```{math}
  :label: loglogfit
  t = Cn^p \qquad \Longrightarrow \qquad \log t = p(\log n) + (\log C).
```

Hence, a graph of $\log t$ as a function of $\log n$ will be a straight line of slope $p$.

## Flops for solving linear equations

Recall how we have proposed to solve the system $\mathbf{A}\mathbf{x}=\mathbf{b}$:

1. Factor $\mathbf{L}\mathbf{U}=\mathbf{A}$ using Gaussian elimination.
1. Solve $\mathbf{L}\mathbf{z}=\mathbf{b}$ for $\mathbf{z}$ using forward substitution.
1. Solve $\mathbf{U}\mathbf{x}=\mathbf{z}$ for $\mathbf{x}$ using backward substitution.

The second and third steps are solved by {ref}`function-forwardsub` and {ref}`function-backsub`. Only one line in each of these functions performs any arithmetic. Take `forwardsub`, for instance. It has a single flop in line 11.  Line 13 computes

```julia  
sum( L[i,j]*x[j] for j=1:i-1 )
```

This line requires $i-1$ multiplications and $(i-2)$ additions, for a total of $2i-3$ flops. Line 14 adds two more flops. These lines are performed within a loop as $i$ ranges from 1 to $n$, so the total count is

```{math}
  :label: trisolveflops
  1 + \sum_{i=1}^n (2i-3) = 1 - 3n + 2 \sum_{i=1}^n.
```

It is not hard to find an exact formula for the sum, but we will use asymptotics to simplify it. Using calculus it can be proved that

```{math}
:label: sumflops
\begin{split}
  \sum_{k=1}^n k&\sim \frac{n^2}{2} = O(n^2), \text{ as $n\to\infty$}, \\
  \sum_{k=1}^n k^2 &\sim \frac{n^3}{3} = O(n^3), \text{ as $n\to\infty$}, \\
  &\vdots \\
  \sum_{k=1}^n k^p &\sim \frac{n^{p+1}}{p+1} = O(n^{p+1}), \text{ as $n\to\infty$},
\end{split}
```

```{margin}
Solving a triangular $n\times n$ system takes $\sim n^2$ flops asymptotically.
```

which holds for any positive integer $p$. This formula has memorable similarity to the antiderivative "power rule" for $\int x^p \, dx$. Applying it to {eq}`trisolveflops` leads us to conclude that solving a triangular linear system of size $n\times n$ takes $\sim n^2$ flops. An analysis of backward substitution yields the same result.

Now let's count flops for LU factorization from the listing of {ref}`function-lufact`. Line 15 requires one division.  Line 16 is an operation on vectors of length $n-j+1$, requiring $n-j+1$ scalar multiplications and the same number of subtractions. These statements occur inside two loops, which requires a nested sum:

```{math}
:label: gecount1
\sum_{j=1}^{n-1} \sum_{i=j+1}^n \bigl[ 2(n-j) +3 \bigr] = \sum_{j=1}^{n-1} (n-j) \bigl[ 2(n-j) +3 \bigr].
```

If we transform the index using $k=n-j$, this becomes

```{margin}
LU factorization of an $n\times n$ matrix takes $\sim\frac{2}{3}n^3$ flops as $n\to \infty$. This dominates the flops for solving an $n\times n$ linear system.
```

```{math}
\begin{split}
  \sum_{k=1}^{n-1} k ( 2k +3 ) &= 2 \sum_{k=1}^{n-1} k^2 + 3 \sum_{k=1}^{n-1} k \\
  &\sim \frac{2}{3}(n-1)^3 + \frac{3}{2}(n-1)^2 \\
  & \sim \frac{2}{3}n^3.
\end{split}
```

```{prf:example} Julia demo
:class: demo
{doc}`demos/flops-lufact`
```

In conclusion, LU factorization takes $\sim\frac{2}{3}n^3$ flops as $n\to \infty$. This dwarfs the $O(n^2)$ count of the triangular system solves.

In practice, flops are not the only aspect of an implementation that occupies significant time. Our position is that counting flops as a measure of performance is a useful oversimplification. We will assume that LU factorization (and as a result, the solution of a linear system of $n$ equations) requires a real-world time that is roughly $O(n^3)$. This growth rate is a great deal more tolerable than, say, $O(2^n)$, but it does mean that for (at this writing) $n$ greater than 10,000 or so, something other than general LU factorization will have to be used.

## Exercises

1. ✍ The following are asymptotic assertions about the limit $n\rightarrow\infty$. In each case, prove the statement true or false.

    **(a)** $n^2 = O(\log n),\quad$ 
    **(b)** $n^{a} = O(n^b)$ if $a\le b,\quad$
    **(c)** $e^n \sim e^{2n},\quad$
    **(d)** $n+\sqrt{n}\sim n+2\sqrt{n}$.

2. ✍ The following are asymptotic assertions about the limit $h\to 0$. In each case, prove the statement true or false.

    **(a)** $h^2\log(h) = O(h^3),\quad$
    **(b)** $h^{a} = O(h^b)$ if $a < b,\quad$
    **(c)** $\sin(h) \sim h,\quad$
    **(d)** $(e^{2h}-1)\sim h$.

    %%
    %% Part (a)
    % False, since $h^2\log(h)/h^3$ in unbounded.
    %% part(b)
    % Since $h^{a-b}$ is unbounded as $h\to 0$ if $a<b$, this is false.
    %% part (c)
    % True, since $\sin(h)/h\to 1$ as $h\to 0$.
    %% part (d)
    % False, since $(e^{2h}-1)/h\to 2$ as $h\to 0$.

3. ✍ Show that the inner product of two $n$-vectors takes exactly $2n-1$ flops.

4. ✍ Show that the multiplication of two $n\times n$ matrices takes $\sim 2n^3$ flops.

5. ✍ This problem is about evaluation of a polynomial $c_1 + c_2 x + \cdots + c_{n}x^{n-1}$.
  
    **(a)** Here is a little code to do the evaluation.

    ``` julia
    y = c[1]
    xpow = 1
    for i = 2:n
        xpow *= x
        y += c[i]*xpow
    end
    ```

    Assuming that `x` is a scalar, how many flops does this function take, as a function of $n$?

    **(b)** Compare the count from (a) to the flop count for Horner's rule, [`horner`](function-horner).
  
    %%
    %% part (a)
    % Line 4 takes 1 flop and line 5 takes 2. So there are $3(n-1)\sim 3n$
    % flops.
    %% part (b)
    % Horner's rule takes $2(n-1)\sim 2n$ flops, or one-third fewer.

6. The exact sums for $p=1,2$ in {eq}`sumflops` are as follows:
  
    ```{math}
    \sum_{k=1}^{n} k = \frac{n(n+1)}{2}, \qquad 
    \sum_{k=1}^{n} k^2 = \frac{n(n+1)(2n+1)}{6}.
    ```

    **(a)** ✍  Use these to find the exact result for {eq}`gecount1`.

    **(b)** ⌨ On one log-log graph, plot the exact expression from part (a) together with the asymptotic result $2n^3/3$ for all $n=10^{1+3i}$, $i=0,\dots,100$.

    **(c)** ⌨ Plot the ratio of the two expressions as a function of $n$ on a plot with a log scale on the $x$-axis and a linear scale on the $y$-axis.
  
    %%
    %% part(a)
    % If you use $k=n-j$, then (2.5.4) becomes $\sum_{k=1}^{n-1} (2k^2+3k)$. Putting in the summation formulas, we get 
    %%
    % $$\frac{3(n-1)n}{2} + \frac{2(n-1)(n)(2n-1)}{6}$$
    %%
    % which becomes 
    %%
    % $$ \frac{1}{6}(9n^2-9n+4n^3-6n^2+2n) = \frac{1}{6}(4n^3+3n^2-7n) $$
    %%
    %% part(b)
   % n = logspace(1,4);
   % exact = (4*n.^3 + 3*n.^2 -7*n)/6;
   % asymp = 2*n.^3/3;
   % loglog(n,[exact(:) asymp(:)])
    %% part (c)
   % semilogx(n,exact./asymp)

7. ✍ Show that for any nonnegative constant integer $m$,
  
    ```{math}
    \sum_{k=0}^{n-m} k^p \sim \frac{n^{p+1}}{p+1}.
    ```
