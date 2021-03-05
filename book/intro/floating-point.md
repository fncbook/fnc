# Floating-point numbers

The real number set $\real$ is infinite in two ways: it is unbounded and continuous. In most practical computing, the second kind of infiniteness is much more consequential than the first kind, so we turn our attention there first. We replace $\real$ with the set $\float$ of {term}`floating point numbers` {index}`floating point numbers`, whose members are zero and all numbers of the form

```{math}
  :label: floatpoint
  \pm (1 + f) \times 2^e,
```

```{index} see: significand; mantissa
```

where $e$ is an integer called the **exponent**, and $1+f$ is the {index}`!mantissa` **mantissa** or **significand**, in which

```{math}
  :label: mantissa
  f = \sum_{i=1}^d b_i \, 2^{-i}, \qquad b_i\in\{0,1\},
```

for a fixed integer $d$. Equation {eq}`mantissa` represents the mantissa as a number in $[1,2)$ in base-2 form. Equivalently,

```{math}
   f = 2^{-d}\, \sum_{i=1}^{d} b_{i} \, 2^{d-i} = 2^{-d} z,
```

for an integer $z$ in the set $\{0,1,\ldots,2^d-1\}$. Consequently, starting at $2^e$ and ending just before $2^{e+1}$ there are exactly $2^d$ evenly spaced numbers belonging to $\float$.

````{prf:example}
  Suppose $d=2$. Hence with $e=0$ in {eq}`floatpoint`, we see that $1$, $5/4$, $3/2$, and $7/4$ are floating point numbers. They are the only members of $\float$ in the half-interval $[1,2)$. Taking $e=1$ gives the floating point numbers in $[2,4)$, specifically $2$, $2.5$, $3$, and $3.5$. Generally the spacing between the elements of $\float$ in $[2^e,2^{e+1})$ is $2^{e-d}$.
````

```{index} see: unit roundoff; machine epsilon
```

Observe that the smallest element of $\float$ that is greater than 1 is $1+2^{-d}$. We define {term}`machine epsilon` (or **machine precision**) as $\macheps = 2^{-d}$.[^macheps]

[^macheps]: The terms machine epsilon, machine precision, and unit roundoff aren't used consistently across references, but the differences are minor for our purposes.

```{margin}
The distance between the floating-point numbers in $\pm[2^e,2^{e+1})$ is $2^e\macheps$.
```

We also suppose the existence of a rounding function $\fl(x)$ that maps every real $x$ to the nearest member of $\float$. If $x$ is positive, we know that it lies in some interval $[2^e,2^{e+1})$, where the spacing between elements of $\float$ is $2^{e-d}$. Therefore we conclude that $|\fl(x)-x| \le \tfrac{1}{2}(2^{e-d})$, which leads to the bound

```{math}
:label: fpbound
\frac{|\fl(x)-x|}{|x|} \le \frac{2^{e-d-1}}{2^e} \le  \tfrac{1}{2}\macheps.
```

Equation {eq}`fpbound` holds true for negative $x$ as well. A mathematically equivalent (see [this exercise](problem-fprelative)) and sometimes more convenient statement is that

```{math}
  :label: fpboundalt
  \fl(x)=x(1+\epsilon), \quad \text{for some $|\epsilon|\le \tfrac{1}{2}\macheps$.}
```

The value of $\epsilon$ depends on $x$, but this dependence is not usually shown explicitly.

## Precision and accuracy

The upshot of floating point representation, as stated in {eq}`fpbound`, is that every real number is represented with a uniformly bounded relative precision. Except for the use of base 2 rather than base 10, floating point representation is a form of scientific notation. For example, Planck's constant is $6.626068\times 10^{-34}$ $\text{m}^2\text{kg}/\text{s}$. If we alter just the last digit from $8$ to $9$, the relative change is

\begin{equation*}
\frac{0.000001\times 10^{-34}}{6.626068\times 10^{-34}} \approx 1.51\times 10^{-7}.
\end{equation*}

This observation justifies the statement that the constant is given with 7 {index}`significant digits` significant digits of precision in base 10. That's in contrast to saying that the value is given to 40 decimal places. A major advantage of floating point is that the relative precision does not depend on the choice of physical units. For instance, when expressed in $\text{eV}\cdot\text{sec}$, Planck's constant is $4.135668\times 10^{-15}$, which has still 7 digits but only 21 decimal places.

It can be easy to confuse the terms **precision** and **accuracy**, especially when looking at the result of a calculation on the computer. The precision of a floating point number is always $d$ binary digits, but not all of those digits may accurately represent an intended value.

Suppose $x$ is a number of interest and $\tilde{x}$ is an approximation to it. The **absolute accuracy** of $\tilde{x}$ is

```{math}
|\tilde{x} - x|,
```

while the **relative accuracy** is

```{math}
\frac{|\tilde{x} - x|}{|x|},
```

Absolute accuracy has the same units as $x$, while relative accuracy is dimensionless. We can also express the relative accuracy as

```{math}
:label: sigdig
\text{accurate digits} = -\log_{10} \left| \frac{\tilde{x}-x}{x} \right|.
```

```{prf:example} Julia demo
:class: demo
{doc}`demos/float-accuracy`
```

We often round this down to an integer, but it does make sense to speak of "almost seven digits" or "ten and a half digits."

## Double precision

```{index} machine epsilon; in double precision
```

```{index} IEEE 754
```

Most numerical computing today is done in the **IEEE 754** standard. This defines **single precision** with $d=23$ binary digits for the fractional part $f$, and the more commonly used {term}`double precision` with $d=52$. In double precision,

```{math}
:label: doubleprec
\macheps = 2^{-52} \approx 2.2\times 10^{-16}.
```

We often speak of double-precision floating point numbers as having about 16 decimal digits. The 52-bit mantissa is paired with a sign bit and 11 binary bits to represent the exponent $e$ in {eq}`floatpoint`, for a total of 64 binary bits per floating point number.

```{prf:example} Julia demo
:class: demo
{doc}`demos/float-julia`
```

Our theoretical description of $\float$ did not place limits on the exponent, but in double precision its range is limited to $-1022\le e \le 1023$. Thus, the largest number is just short of $2^{1024}\approx 2\times 10^{308}$, which is more than enough in most applications. Results that should be larger are said to {index}`overflow` and will actually result in the value `Inf`. Similarly, the smallest positive number is $2^{-1022}\approx 2\times 10^{-308}$, and smaller values are said to {index}`underflow` to zero.[^denormalized]

[^denormalized]: Actually, there are some still-smaller **denormalized** numbers that have less precision, but we won't use that level of detail.

Note the crucial difference between $\macheps=2^{-52}$, which is the distance between 1 and the next larger double-precision number, and $2^{-1022}$, which is the smallest positive double-precision number. The former has to do with relative precision, while the latter is about absolute precision. Getting close to zero always requires a shift in thinking to absolute precision, because any finite error is infinite relative to zero.

One more double precision value is worthy of note: {term}`NaN`, which stands for **Not a Number**. It is the result of an undefined arithmetic operation such as 0/0.

## Floating point arithmetic

```{index} floating-point numbers
```

Computer arithmetic is performed on floating-point numbers and returns floating-point results. We assume the existence of machine-analog operations for real functions such as $+$, $-$, $\times$, $/$, $\sqrt{\quad}$, and so on. Without getting into the details, we will suppose that each elementary machine operation creates a floating point result whose relative error is bounded by $\macheps$. For example, if $x$ and $y$ are in $\float$, then for machine addition $\oplus$ we have the bound

```{math}
\frac{|(x \oplus y)-(x+y)|}{|x+y|} \le \macheps.
```

```{prf:example} Julia demo
:class: demo
{doc}`demos/float-arithmetic`
```

Hence the relative error in arithmetic is practically the same as for the floating point representation itself. However, even playing by these rules can lead to disturbing results.

There are two ways to look at the result in {doc}`demos/float-arithmetic`. On one hand, its two versions of the result differ by less than $1.2\times 10^{-16}$, which is very small—not just in everyday terms, but with respect to the operands, which are all close to 1 in absolute value. On the other hand, the difference is as large as the exact result itself! We formalize and generalize this observation in the next section. In the meantime, keep in mind that exactness cannot be taken for granted in floating point computation. Even ideally, we should not expect that two mathematically equivalent results will be equal, only that they be relatively close together.

## Exercises

```{note}
Exercises marked with ✍ are intended to be done by hand or with the aid of a simple calculator. Exercises marked with ⌨ are intended to be solved by using a computer.
```

1. ✍ Consider a floating point set $\float$ defined by {eq}`floatpoint` and {eq}`mantissa` with $d=4$.
  
    **(a)** How many elements of $\float$ are there in the real interval $[1/2,4]$, including the endpoints?

    **(b)** What is the element of $\float$ closest to the real number $1/10$?

    **(c)** What is the smallest positive integer not in $\float$?

    (problem-fprelative)=

2. ✍ Prove that {eq}`fpbound` is equivalent to {eq}`fpboundalt`. (This means showing first that {eq}`fpbound` implies {eq}`fpboundalt`, and then separately that {eq}`fpboundalt` implies {eq}`fpbound`.)

3. ⌨ There are much better rational approximations to $\pi$ than $22/7$. For each one below, find its absolute and relative accuracy, and (rounding down to an integer) the number of accurate digits. (In Julia, the variable `pi` is set to a floating-point approximation to $\pi$ by default.)

    **(a)** 355/113

    **(b)** 103638/32989

    ```{index} IEEE 754
    ```

4. ✍ IEEE 754 **single precision** specifies that 23 binary bits are used for the mantissa $f$ in {eq}`mantissa`. Because they need less storage space and can be operated on more quickly than double precision values, single precision values can be useful in low-precision applications.

    **(a)** In base-10 terms, what is the first single precision number greater than $1$ in this system?

    **(b)** What is the smallest positive integer that is not a single precision number?

5. ⌨ Julia defines a function `nextfloat` that gives the next-larger floating-point value of a given number. What is the next float past `floatmax`?

6. ⌨ (This problem requires the use of loops in Julia.) It's reasonable to expect that floating point errors accumulate randomly during a long computation, creating what is known as a **random walk**. On average we expect as many errors to be negative as positive, so they tend to partially cancel out. Suppose we define a random sequence by $x_0=0$ and $x_{n}=x_{n-1}\pm 1$ for $n\ge 1$, with the signs chosen by tossing a fair coin for each $n$. Let $\alpha_n$ and $\beta_n$ be the average value of $x_n$ and $|x_n|$, respectively, over all such walks. Then a classic result of probability is that $\alpha_n=0$ and

    ```{math}
    \lim_{n\to\infty} \frac{\pi \beta_n^2}{2n}=1.
    ```

    In Julia the function `randn` simulates drawing numbers from the normal or Gaussian distribution (i.e., the bell curve) with mean zero and variance 1. Choose a unique positive integer seed value $s$ (for example, use the last 5 digits of your phone number) and enter `import Random; Random.seed!(s)` to initialize the random number generator. Then the following code generates one random walk for $n=10^4$:

    ``` julia
    r = randn(10000)           # draw random numbers
    x = sum(sign(z) for z in r)
    ```

    Perform a million random walks, computing the average values of $x_{10000}$ and $|x_{10000}|$. Compare these to $\alpha_n=0$ and $\beta_n\approx \sqrt{2n/\pi}$ at $n=10000$.
