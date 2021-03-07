# Zero-stability of multistep methods

```{index} multistep method
```

````{prf:example} Julia demo
:class: demo
:label: demos-zs-LIAF
{doc}`demos/zs-LIAF`
````

For one-step methods such as Runge--Kutta, {prf:ref}`theorem-onestepGTE` guarantees that the method converges and that the global error is of the same order as the local truncation error. For multistep methods, however, a new wrinkle is introduced. As an example, it can be checked that the 2-step method "LIAF," defined by

```{math}
  :label: LIAF
  \mathbf{u}_{i+1} = -4u_i + 5u_{i-1} + h(4f_i + 2f_{i-1}),
```

is third-order accurate. Yet it is not useful, as can be seen in {prf:ref}`demos-zs-LIAF`.

```{index} instability; of a mutistep formula
```

The source of the instability in {prf:ref}`demos-zs-LIAF` is not hard to identify. First though, we're going to encounter the possibility of complex numbers shortly. Because it's easy to confuse the node index $i$ with the imaginary unit, we'll switch to using $m$ as the step index in this section.

Let's recall that we can rewrite {eq}`LIAF` as $\rho(\mathcal{Z})u_{m-1}=h \sigma(\mathcal{Z})u_{m-1}$ using the forward shift operator $\mathcal{Z}$:

```{math}
  :label: LIAFshift
  (\mathcal{Z}^2 + 4\mathcal{Z} - 5) u_{m-1} = h(4\mathcal{Z} + 2)f_{m-1}.
```

(See {eq}`multistepshift`.) Next, suppose that $h$ is negligible in {eq}`LIAFshift`. Then the numerical solution of LIAF is roughly defined by

```{math}
  :label: LIAFshiftzero
  (\mathcal{Z}^2 + 4\mathcal{Z} - 5) u_{m-1} = 0.
```

The graph in {prf:ref}`demos-zs-LIAF` strongly suggests that for small $h$, $|u_m|\approx c \alpha^m$ for some $\alpha>1$ as $m$ gets large. So we are motivated to posit $u_m = c z^m$ for all $m$ and see if we can prove that it is an exact solution. The beauty of this choice is that $\mathcal{Z} u_i = z u_i$; that is, the "shift ahead" operator on the sequence $u_0,u_1,\ldots$ is identical to "multiply by $z$." This observation implies that we can solve {eq}`LIAFshiftzero` with $u_m = c z^m$ for all $m$ simply by finding a numerical value of $z$ such that

```{math}
  :label: LIAFcharzero
  z^2 + 4z - 5 = 0.
```

We've arrived at a polynomial rootfinding problem! The roots for this method are $z=1$ and $z=-5$.

Now we see why exponentially large numbers were observed in {prf:ref}`demos-zs-LIAF`: the sequence $u_m=c (-5)^m$ is approximately a solution of the multistep formula for small values of $h$. This causes exponential growth that drowns out the solution we were trying to find. In the general case, whenever there is a root $r$ of $\rho(z)=0$ such that $|r|>1$, we should expect exponentially growing solutions as $h\to 0$ and $m\to\infty$.

## The root condition

```{margin}
Zero-stability requires that as $h\to 0$, every numerical solution produced by the multistep formula remains bounded throughout the time interval.
```

```{index} zero-stability
```

The property we lacked in {prf:ref}`demos-zs-LIAF` is called **zero-stability**. To state it precisely, zero-stability requires that as $h\to 0$, every numerical solution produced by the multistep formula remains bounded throughout $a\le t_m \le b$. Without this property, any kind of error, whether from truncation or roundoff, will get exponentially amplified and overwhelm convergence to the exact solution. The following theorem concisely summarizes when we can expect zero-stability.

````{prf:criterion} Root condition
:label: theorem-rootcondition
A linear multistep method is zero-stable (i. e., has only bounded solutions as $h\rightarrow 0$) if and only if every root $r$ of the generating polynomial $\rho(z)$ satisfies $|r|\le 1$, and any root $r$ with $|r|=1$ is simple (that is, $\rho'(r)\neq 0$).
````

````{prf:proof}
(Partial proof) As described above, the values produced by the numerical method approach solutions of the difference equation $\rho(\mathcal{Z})u_m=0$. We consider only the case where the roots $r_1,\ldots,r_k$ of $\rho(z)$ are all distinct. Then $u_m=(r_j)^m$ is a solution of $\rho(\mathcal{Z})u_m=0$ for each $j=1,\ldots,k$. By linearity,

```{math}
  u_m = c_1 (r_1)^m + c_2 (r_2)^m + \cdots + c_k (r_k)^m
```

is a solution for any values of $c_1,\ldots,c_k$. These constants are determined uniquely by the starting values $u_0,\ldots,u_{k-1}$ (we omit the proof). Now, if all the roots satisfy $|r_j|\le 1$, then

```{math}
  |u_m| \le \sum_{j=1}^k |c_j| |r_j|^m \le \sum_{j=1}^k |c_j|,
```

independently of $h$ and $m$. This proves zero-stability. Conversely, if some $|r_j|>1$, then $|u_m|$ cannot be bounded above by a constant independent of $m$. Since $b=t_m$, $m\to\infty$ at $t=b$ as $h\to 0$, so zero-stability cannot hold.
````

````{prf:example}
A $k$-step Adams method has $\rho(z) = z^k - z^{k-1} = z^{k-1}(z-1)$. Hence 1 is a simple root and 0 is a root of multiplicity $k-1$. So the Adams methods are all stable.
````

## Dahlquist theorems

It turns out that lacking zero-stability is the only thing that can go wrong for a consistent multistep method. (Recall that a method is consistent if its local truncation error is $O(h)$ as $h\to 0$.)

```{prf:theorem} Dahlquist equivalence
:label: theorem-dahlequiv
A linear multistep method is convergent if and only if it is consistent and zero-stable.
```

The Dahlquist equivalence theorem is one of the most important and celebrated in the history of numerical analysis. It can be proved more precisely that a zero-stable, consistent method is convergent in the same sense as {prf:ref}`theorem-onestepGTE`, with the error between numerical and exact solutions being of the same order as the local truncation error, for a wide class of problems.

You may have noticed that the Adams and BD formulas use only about half of the available data from the past $k$ steps, i.e., they have many possible coefficients set to zero. For instance, a $k$-step AB method uses only the $f_j$-values and has order $k$. The order could be made higher by also using $u_j$-values, like the LIAF method does for $k=2$. Also like the LIAF method, however, such attempts are doomed by instability.

````{prf:theorem} First Dahlquist stability barrier
The order of accuracy $p$ of a stable $k$-step linear multistep method satisfies

```{math}
p \le
\begin{cases}
  k+2 & \text{if $k$ is even},\\
  k+1 & \text{if $k$ is odd},\\
  k & \text{if the method is explicit.}
\end{cases}
```
````

## Exercises

1. ✍ Show that the LIAF method {eq}`LIAF` has order of accuracy equal to three.

    ````{only} solutions
    ````

2. ✍ / ⌨  Verify that the order of accuracy of the given multistep method is at least one. Then apply {prf:ref}`theorem-rootcondition` to determine whether it is zero-stable.

    **(a)** BD2

    **(b)** BD3

    **(c)** $u_{i+1}=u_{i-1}+2hf_i$

    **(d)** $u_{i+1} = -u_i +u_{i-1} + u_{i-2} + \frac{2h}{3}(4f_i+f_{i-1}+f_{i-2})$

    **(e)** $u_{i+1} = u_{i-3} + \frac{4h}{3} ( 2f_i - f_{i-1} + 2f_{i-2})$

    **(f)** $u_{i+1} = -2u_i + 3u_{i-1} + h (f_{i+1}+2f_i+f_{i-1})$
  
    ````{only} solutions
    ````

3. ✍  A Fibonacci sequence is defined by $u_{i+1}=u_i+u_{i-1}$, where $u_0$ and $u_1$ are seed values. Using the proof of {prf:ref}`theorem-rootcondition`, find $r_1$ and $r_2$ such that $u_i=c_1(r_1)^i+c_2(r_2)^i$ for all $i$.

4. ✍ Suppose that $r$ is a root of multiplicity at least two for $\rho(z)$. Show that $u_m = c m r^m$ is a solution of the difference equation $\rho(\mathcal{Z})u_m=0$. 

    ````{only} solutions
    ````
