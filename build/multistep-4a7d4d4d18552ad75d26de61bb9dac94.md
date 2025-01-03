---
numbering:
  enumerator: 6.6.%s
---
(section-ivp-multistep)=
# Multistep methods

In Runge–Kutta methods we start at $u_i$ to find ${u}_{i+1}$, taking multiple $f$-evaluations (stages) to achieve high accuracy. In contrast, multistep methods boost accuracy by employing more of the history of the solution, taking information from the recent past. For the discussion in this and following sections, we introduce the shorthand notation

```{math}
:label: multistepfi
f_i = f(t_i,u_i).
```

```{index} ! multistep method, ! implicit IVP solver
```

::::{prf:definition} Multistep method for IVPs
A $k$-step **multistep** (or *linear multistep*) method is given by the difference equation

```{math}
:label: multistep
\begin{split}
u_{i+1} &= a_{k-1}u_i + \cdots + a_0 u_{i-k+1} \qquad \\ 
& \qquad + h ( b_kf_{i+1} + \cdots + b_0 f_{i-k+1}),
\end{split}
```

where the $a_j$ and the $b_j$ are constants. If $b_k=0$, the method is **explicit**; otherwise, it is **implicit**.  
::::

The quantities $u$ and $f$ in {eq}`multistep` are shown as scalars, but in general they can be vectors. 

In order to use {eq}`multistep` as a numerical method, we iterate through $i=k-1,\ldots,n-1$. The value $u_0$ is determined by the initial condition, but we also need some way of generating the **starting values**

```{math}
u_1=\alpha_1, \quad \ldots \quad u_{k-1}=\alpha_{k-1}.
```

```{index} Runge–Kutta method
```

In practice the starting values are often found using an RK formula.[^whyRK]

[^whyRK]: If we must use an RK method to start anyway, why bother with multistep formulas at all? The answer is that multistep methods can be more efficient in some problems, even at the same order of accuracy.

The difference formula {eq}`multistep` defines ${u}_{i+1}$ in terms of known values of the solution and its derivative from the past. In the explicit case with $b_k=0$, Equation {eq}`multistep` immediately gives a formula for the unknown quantity ${u}_{i+1}$ in terms of values at time level $t_i$ and earlier. Thus only one new evaluation of $f$ is needed to make a time step, provided that we store the recent history. 

For an implicit method, however, $b_k\neq 0$ and {eq}`multistep` has the form

```{math}
  {u}_{i+1} - hb_kf(t_{i+1},{u}_{i+1}) = F(u_i,u_{i-1},\ldots,u_{i-k+1}).
```

```{index} rootfinding problem
```

Now the unknown ${u}_{i+1}$ that we seek appears inside the function $f$. In general this equation is a nonlinear rootfinding problem for ${u}_{i+1}$ and is not solvable in a finite number of steps by a formula. The implementation of both explicit and implicit multistep formulas is discussed in detail in {numref}`section-ivp-implicit`.


```{index} ! Adams--Bashforth formula, ! Adams--Moulton formula, ! backward differentiation formula for IVPs
```

As with RK formulas, a multistep method is entirely specified by the values of a few constants. {numref}`table-adams` and {numref}`table-BD` present some of the most well-known and important formulas. The **Adams–Bashforth** (AB) methods are explicit, while **Adams–Moulton** (AM) and **backward differentiation formulas** (BD) are implicit. The tables also list the methods' order of accuracy, to be defined shortly. We adopt the convention of referring to a multistep method by appending its order of accuracy to a two-letter name abbreviation, e.g., the AB3 method.

```{index} ! trapezoid formula; for an IVP
```

(table-adams)=

```{table} Coefficients of Adams multistep formulas. All have $a_{k-1}=1$ and $a_{k-2} = \cdots = a_0 = 0$.

| name/order | steps $k$ |       $b_k$       |     $b_{k-1}$     |     $b_{k-2}$      |     $b_{k-3}$     |     $b_{k-4}$     |
|:----------:|:---------:|:-----------------:|:-----------------:|:------------------:|:-----------------:|:-----------------:|
|    AB1     |     1     |         0         |         1         |      (Euler)       |                   |                   |
|    AB2     |     2     |         0         |   $\frac{3}{2}$   |   $-\frac{1}{2}$   |                   |                   |
|    AB3     |     3     |         0         |  $\frac{23}{12}$  |  $-\frac{16}{12}$  |  $\frac{5}{12}$   |                   |
|    AB4     |     4     |         0         |  $\frac{55}{24}$  |  $-\frac{59}{24}$  |  $\frac{37}{24}$  |  $-\frac{9}{24}$  |
|    AM1     |     1     |         1         | (Backward Euler)  |                    |                   |                   |
|    AM2     |     1     |   $\frac{1}{2}$   |   $\frac{1}{2}$   |    (Trapezoid)     |                   |
|    AM3     |     2     |  $\frac{5}{12}$   |  $\frac{8}{12}$   |  $-\frac{1}{12}$   |                   |                   |
|    AM4     |     3     |  $\frac{9}{24}$   |  $\frac{19}{24}$  |  $-\frac{5}{24}$   |  $\frac{1}{24}$   |                   |
|    AM5     |     4     | $\frac{251}{720}$ | $\frac{646}{720}$ | $-\frac{264}{720}$ | $\frac{106}{720}$ | $-\frac{19}{720}$ |
```

(table-BD)=

```{table} Coefficients of backward differentiation formulas. All  have $b_k\neq 0$ and $b_{k-1} = \cdots = b_0 = 0$.

| name/order | steps $k$ | $a_{k-1}$       | $a_{k-2}$        | $a_{k-3}$       | $a_{k-4}$       | $b_k$           |
|:----------:|:---------:|:---------------:|:----------------:|:---------------:|:---------------:|:---------------:|
| BD1        | 1         | 1               | (Backward Euler) |                 |                 | 1               |
| BD2        | 2         | $\frac{4}{3}$   | $-\frac{1}{3}$   |                 |                 | $\frac{2}{3}$   |
| BD3        | 3         | $\frac{18}{11}$ | $-\frac{9}{11}$  | $\frac{2}{11}$  |                 | $\frac{6}{11}$  |
| BD4        | 4         | $\frac{48}{25}$ | $-\frac{36}{25}$ | $\frac{16}{25}$ | $-\frac{3}{25}$ | $\frac{12}{25}$ |
```

## Generating polynomials

```{index} ! generating polynomials
```

An alternative description of a multistep method is the **generating polynomials**

```{math}
:label: genpoly
\begin{split}
  \rho(z) &= z^k - a_{k-1} z^{k-1} - \cdots - a_0,\\
  \sigma(z) &= b_k z^k + b_{k-1}z^{k-1} + \cdots + b_0.
\end{split}
```

For example, the AB3 method is completely specified by

```{math}
  \rho(z) = z^3-z^2, \qquad \sigma(z) = \tfrac{1}{12}(23z^2-16z+5).
```

:::{prf:observation}
Let $\rho$ and $\sigma$ be the generating polynomials of a multistep method. Then:

1. The polynomial $\rho(z)$ is monic (i.e., its leading term has a unit coefficient).
2. The degree of $\rho$ is the number of steps $k$. 
3. The degree of $\sigma(z)$ is $k$ for an implicit method and less than $k$ for an explicit method. 
:::

The connection between the generating polynomials and the numerical method requires a little abstraction. Let $\mathcal{Z}$ be a *forward-shift operator*, so that, for example, $\mathcal{Z} t_i = t_{i+1}$, $\mathcal{Z}^3 u_{i-1} = u_{i+2}$, and so on. With this, the difference formula {eq}`multistep` can be written concisely as

```{math}
:label: multistepshift
  \rho(\mathcal{Z}) u_{i-k+1} = h \sigma(\mathcal{Z}) f_{i-k+1}.
```

## Truncation and global error

```{index} ! truncation error; of a multistep IVP formula
```

The definition of local truncation error is easily extended to multistep methods. 

(definition-multistep-lte)=
:::{prf:definition} LTE and order of accuracy for a multistep IVP method
For the multistep formula defined by {eq}`multistep`, the **local truncation error** is 

```{math}
:label: MSLTE
  \tau_{i+1}(h) = \frac{\hat{u}(t_{i+1}) - a_{k-1}\hat{u}(t_i) - \cdots - a_0
    \hat{u}(t_{i-k+1})}{h} - \bigl[ b_kf(t_{i+1},\hat{u}(t_{i+1})) + \cdots +
  b_0f(t_{i-k+1},\hat{u}(t_{i-k+1})) \bigr].
```
If the local truncation error satisfies $\tau_{i+1}(h)=O(h^p)$, then $p$ is the **order of accuracy** of the formula. If $p>0$, the method is **consistent**.
::::

````{prf:example}
The first-order Adams–Moulton method is also known as **backward Euler**, because its difference equation is

```{math}
  {u}_{i+1} = u_i + hf_{i+1},
```

which is equivalent to a backward-difference approximation to $u'(t_{i+1})$. AM1 is characterized by $\rho(z) = z-1$ and $\sigma(z) = z$.

To derive the LTE, we use the definition:

```{math}
\begin{split}
    h\tau_{i+1}(h) &= \hat{u}(t_{i+1}) - \hat{u}(t_i) - hf\bigl(t_{i+1},\hat{u}(t_{i+1})\bigr) \\
    &= \hat{u}(t_i) + h\hat{u}'(t_i) + \frac{h^2}{2}\hat{u}''(t_i) + O(h^3)
    - \hat{u}(t_i) -h \hat{u}'(t_{i+1}) \\
    &= h\hat{u}'(t_i) + \frac{h^2}{2}\hat{u}''(t_i) + O(h^3)
    - h[\hat{u}'(t_i) + h\hat{u}''(t_i) + O(h^2)]\\
    &= - \frac{h^2}{2}\hat{u}''(t_i) + O(h^3).
\end{split}
```

Thus $\tau_{i+1}(h)=O(h)$ and AM1 (backward Euler) is a first-order method.
````

````{prf:example}
The AB2 method has the formula

```{math}
  {u}_{i+1} = u_i + h\left(\frac{3}{2} f_i - \frac{1}{2} f_{i-1} \right).
```

The generating polynomials are $\rho(z)=z^2-z$ and $\sigma(z) = (3z-1)/2$. We find that the method is second order from the LTE:

```{math}
\begin{split}
  h\tau_{i+1}(h)
  & = \hat{u}(t_{i+1}) - \hat{u}(t_i) - h\left[
    \frac{3}{2}f(t_i,\hat{u}(t_i)) - \frac{1}{2}f(t_{i-1},\hat{u}(t_{i-1}))
    \right]                                                                                   \\
  & = \hat{u}(t_i) + h\hat{u}'(t_i) + \frac{h^2}{2}\hat{u}''(t_i) + \frac{h^3}{6}\hat{u}'''(t_i) + O(h^4) \\
  & \qquad - \hat{u}(t_i) - \frac{3h}{2}\hat{u}'(t_i)  \\
  &\qquad  + \frac{h}{2} \bigl[\hat{u}'(t_i) - h\hat{u}''(t_i) + \frac{h^2}{2}\hat{u}'''(t_i) + O(h^3)\bigr]        \\
  & = \frac{5h^3}{12}\hat{u}'''(t_i) + O(h^4),
\end{split}
```

so that $\tau_{i+1}(h)=O(h^2)$.
````

Although we will not present the analysis, the main conclusion for the multistep methods in this section is the same as for one-step methods.

:::{prf:observation}
The global error of each method in {numref}`table-adams` and {numref}`table-BD` converges at the same order as the local truncation error.
:::

## Derivation of the formulas

Where do coefficients like those in {numref}`table-adams` come from? There are different ways to answer that question, but Adams and BD methods have distinctive stories to tell. The derivation of Adams methods begins with the observation that

```{math}
:label: adamsderive
  \hat{u}(t_{i+1}) = \hat{u}(t_i) + \int_{t_i}^{t_{i+1}} \hat{u}'(t) \, dt =
  \hat{u}(t_i) + \int_{t_i}^{t_{i+1}} f\bigl(t,\hat{u}(t)\bigr) \, dt.
```

As a result, a $k$-step Adams method always has $\rho(z)=z^k-z^{k-1}$. While the integrand above is unknown over the interval of integration, we can approximate it by a polynomial interpolant of the historical values of $f$. That polynomial can be integrated analytically, leading to a derivation of the coefficients $b_0,\ldots,b_k$.

(example-am2derive)=

````{prf:example}
Let's derive a one-step AM method using the two values $(t_i,f_i)$ and $(t_{i+1},f_{i+1})$. The interpolating polynomial is the linear function

```{math}
p(t) = f_i\frac{t_{i+1}-t}{t_{i+1}-t_i} + f_{i+1}\frac{t-t_i}{t_{i+1}-t_i}.
```

Things become a little easier with the change of variable $s=t-t_i$ and applying $h=t_{i+1}-t_i$:

```{math}
\int_{t_i}^{t_{i+1}} p(t)  \, d t = \int_0^h p(t_i+s) \, d s
= h^{-1} \int_0^h [ (h-s)f_i + s f_{i+1} ]\, d s = \frac{h}{2}(f_i + f_{i+1}).
```

```{index} trapezoid formula; for an IVP
```

Hence $\sigma(z)=\tfrac{1}{2}z + \tfrac{1}{2}$. Like the trapezoid formula for a definite integral, AM2 computes the exact integral of a piecewise linear interpolant, and it often goes by the *trapezoid* name as well.
````

In AB methods, the interpolating polynomial has degree $k-1$, which means that its interpolation error is $O(h^k)$. Upon integrating we get a local error of $O(h^{k+1})$, which reduces to a global error of $O(h^k)$. The AM interpolating polynomial is one degree larger, so its order of accuracy is one higher for the same number of steps.

The idea behind backward differentiation formulas is complementary to that for Adams: Interpolate solution values ${u}_{i+1},\ldots,u_{i-k+1}$ by a polynomial $q$, and then, motivated by $f(t,\hat{u})=\hat{u}'(t)$, set

```{math}
:label: bdfderive
  f_{i+1} =q'(t_{i+1}).
```

The quantity $q'(t_{i+1})$ can be approximated by a finite difference of the past solution values, leading to the coefficients of $\rho(z)$ and $\sigma(z)=b_k z^k$.  

:::{prf:example} 
Consulting {numref}`table-FDforward`, we find the finite-difference approximation

$$
q'(t_{i+1}) \approx \frac{1}{h} \left( \frac{3}{2} u_{i+1} - 2 u_i + \frac{1}{2} u_{i-1} \right),
$$

from which we get

$$
h f_{i+1} = \frac{3}{2} u_{i+1} - 2 u_i + \frac{1}{2} u_{i-1}.
$$

Rearranging and Normalizing by the coefficient of $u_{i+1}$ gives $\rho(z)=z^2 + \tfrac{4}{3}z - \tfrac{1}{3}$ and $\sigma(z) = \tfrac{2}{3}z^2$, which is the BD2 method.
:::

## Exercises

1. ✍ For each method, write out the generating polynomials $\rho(z)$ and $\sigma(z)$.

    **(a)** AM2,
    **(b)** AB2,
    **(c)** BD2,
    **(d)** AM3,
    **(e)** AB3.

2. ✍ Write out by hand an equation that defines the first solution value $u_1$ produced by AM1 (backward Euler) for each IVP. (Reminder: This is an implicit formula.)

    **(a)** $u' = -2t u, \quad 0 \le t \le 2, \quad u_0 = 2, \quad h = 0.2$

    **(b)** $u' = u + t, \quad 0 \le t \le 1, \quad u_0 = 2, \quad h = 0.1$

    **(c)** $(1+x^3)uu' = x^2,\quad 0 \le x \le 3, \quad u_0=1, , \quad h = 0.5$

3. ✍ Do the preceding exercise for AM2 (trapezoid) instead of backward Euler.

4. ✍ For each method, find the leading term in the local truncation error using {eq}`MSLTE`.

    **(a)** AM2,
    **(b)** AB2,
    **(c)** BD2.
  
5. ✍/ ⌨ For each method, find the leading term in the local truncation error using {eq}`MSLTE`. (Computer algebra is recommended.)

    **(a)** AM3,
    **(b)** AB3,
    **(c)** BD4.
  
6. ✍ A formula for the quadratic polynomial interpolant through the points $(s_1,y_1)$, $(s_2,y_2)$, and $(s_3,y_3)$ is
  
    ```{math}
    p(x) = \frac{(x-s_2)(x-s_3)}{(s_1-s_2)(s_1-s_3)}\,y_1 +
            \frac{(x-s_1)(x-s_3)}{(s_2-s_1)(s_2-s_3)}\,y_2 +
            \frac{(x-s_1)(x-s_2)}{(s_3-s_1)(s_3-s_2)}\,y_3.
    ```

    **(a)** Use {eq}`adamsderive` and a polynomial interpolant through three points to derive the coefficients of the AM3 method.

    **(b)** Use {eq}`bdfderive` and a polynomial interpolant through three points to derive the coefficients of the BD2 method.
  
7. ✍ By doing series expansion about the point $z=1$, show for BD2 that
  
    ```{math}
    \frac{\rho(z)}{\sigma(z)} - \log(z-1) = O\bigl( (z-1)^3 \bigr).
    ```

8. ✍/ ⌨  By doing series expansion about the point $z=1$, show for AB3 and AM3 that
  
    ```{math}
    \frac{\rho(z)}{\sigma(z)} - \log(z-1) = O\bigl( (z-1)^4 \bigr).
    ```

    (Computer algebra is recommended.)


