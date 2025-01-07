---
numbering:
  enumerator: 1.1.%s
---
(section-intro-floating-point)=
# Floating-point numbers

The real number set $\real$ is infinite in two ways: it is unbounded and continuous. In most practical computing, the latter kind of infiniteness is much more consequential than the former, so we turn our attention there first. 

```{index} ! floating-point numbers
```

::::{prf:definition} Floating-point numbers
The set $\float$ of **floating-point numbers** consists of zero and all numbers of the form

```{math}
:label: floatpoint
  \pm (1 + f) \times 2^n,
```

```{index} see: mantissa; significand
```
```{index} ! significand
```

where $n$ is an integer called the **exponent**, and $1+f$ is the **significand**, in which

```{math}
:label: mantissa
  f = \sum_{i=1}^d b_i \, 2^{-i}, \qquad b_i\in\{0,1\},
```

for a fixed integer $d$ called the binary **precision**. 
::::

Equation {eq}`mantissa` represents the significand as a number in $[1,2)$ in base-2 form. Equivalently,

```{math}
:label: fp-integer
f = 2^{-d}\, \sum_{i=1}^{d} b_{i} \, 2^{d-i} = 2^{-d} z
```

for an integer $z$ in the set $\{0,1,\ldots,2^d-1\}$. Consequently, starting at $2^n$ and ending just before $2^{n+1}$ there are exactly $2^d$ evenly spaced numbers belonging to $\float$. 

(example-float-2bits)=
````{prf:example}
Suppose $d=2$. Taking $n=0$ in {eq}`floatpoint`, we enumerate

$$
1 + \frac{0}{4}, \: 1 + \frac{1}{4}, \: 1 + \frac{2}{4}, \: 1 + \frac{3}{4}.
$$

These are the only members of $\float$ in the semi-closed interval $[1,2)$, and they are separated by spacing $\tfrac{1}{4}$. 

Taking $n=1$ doubles each of the values in the list above, and $n=-1$ halves them. These give the floating-point numbers in $[2,4)$ and $[1/2,1)$, respectively. The spacing between them also is doubled and halved, respectively. 
````

Observe that the smallest element of $\float$ that is greater than 1 is $1+2^{-d}$, and we call the difference *machine epsilon*.[^macheps]

[^macheps]: The terms machine epsilon, machine precision, and unit roundoff aren't used consistently across references, but the differences are not consequential for our purposes.

```{index} ! machine epsilon, unit roundoff 
```

::::{prf:definition} Machine epsilon
For a floating-point set with $d$ binary digits of precision, **machine epsilon** (or *machine precision*) is $\macheps = 2^{-d}$.
::::

We define the rounding function $\fl(x)$ as the map from real number $x$ to the nearest member of $\float$. The distance between the floating-point numbers in $[2^n,2^{n+1})$ is $2^n\macheps=2^{n-d}$. As a result, every real $x \in [2^n,2^{n+1})$ is no farther than $2^{n-d-1}$ away from a member of $\float$. Therefore we conclude that $|\fl(x)-x| \le \tfrac{1}{2}(2^{n-d})$, which leads to the bound

```{math}
:label: fpbound
\frac{|\fl(x)-x|}{|x|} \le \frac{2^{n-d-1}}{2^n} \le  \tfrac{1}{2}\macheps.
```

In words, every real number is represented with a uniformly bounded relative precision. Inequality {eq}`fpbound` holds true for negative $x$ as well. In [Exercise 2](#problem-fp-fprelative) you are asked to show that an equivalent statement is that

```{math}
:label: fpboundalt
  \fl(x)=x(1+\epsilon) \quad \text{for some $|\epsilon|\le \tfrac{1}{2}\macheps$.}
```

The value of $\epsilon$ depends on $x$, but this dependence is not usually shown explicitly.

## Precision and accuracy

```{index} significant digits
``` 

It may help to recast {eq}`floatpoint` and {eq}`mantissa` in terms of base 10:

$$
\pm \left(b_0 + \sum_{i=1}^d b_i \, 10^{-i} \right) \times 10^n = \pm (b_0.b_1b_2\cdots b_d) \times 10^n,
$$

where each $b_i$ is in $\{0,1,\ldots,9\}$ and $b_0\neq 0$. This is simply scientific notation with $d+1$ significant digits. For example, Planck's constant is $6.626068\times 10^{-34}$ m${}^2\cdot$kg/sec to seven digits. If we alter just the last digit from $8$ to $9$, the relative change is

$$
\frac{0.000001\times 10^{-34}}{6.626068\times 10^{-34}} \approx 1.51\times 10^{-7}.
$$

We therefore say that the constant is given with 7 decimal digits of precision. That's in contrast to noting that the value is given to 40 decimal *places*. A major advantage of floating point is that the relative precision does not depend on the choice of physical units. For instance, when expressed in eV$\cdot$sec, Planck's constant is $4.135668\times 10^{-15}$, which still has 7 digits but only 21 decimal places.

Floating-point precision functions the same way, except that computers prefer base 2 to base 10. The **precision** of a floating-point number is always $d$ binary digits, implying a resolution of the real numbers according to {eq}`fpbound`. 

```{index} ! accuracy (relative vs. absolute)
```

It can be easy to confuse precision with **accuracy**, especially when looking at the result of a calculation on the computer. Every result is computed and represented using $d$ binary digits, but not all of those digits may accurately represent an intended value. Suppose $x$ is a number of interest and $\tilde{x}$ is an approximation to it. The **absolute accuracy** of $\tilde{x}$ is

```{math}
|\tilde{x} - x|,
```

while the **relative accuracy** is

```{math}
\frac{|\tilde{x} - x|}{|x|}.
```

Absolute accuracy has the same units as $x$, while relative accuracy is dimensionless. We can also express the relative accuracy as the **number of accurate digits**, computed in base 10 as

```{math}
:label: sigdig
-\log_{10} \left| \frac{\tilde{x}-x}{x} \right|.
```

We often round this value down to an integer, but it does make sense to speak of "almost seven digits" or "ten and a half digits."

(demo-float-accuracy)= 
:::::::{prf:example} Absolute and relative accuracy
``````{tab-set} 
`````{tab-item} Julia
:sync: julia
::::{dropdown} ""
:open:
:::{include} julia/accuracy.ipynb
:::
::::
`````

`````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-float-accuracy-matlab
:::
````` 

`````{tab-item} Python
:sync: python
:::{embed} #demo-float-accuracy-python
:::
````` 
``````
:::::::


## Double precision

```{index} machine epsilon; in double precision, IEEE 754
```

Most numerical computing today is done in the **IEEE 754** standard. This defines **single precision** with $d=23$ binary digits for the fractional part $f$, and the more commonly used **double precision** with $d=52$. In double precision,

```{math}
:label: doubleprec
\macheps = 2^{-52} \approx 2.2\times 10^{-16}.
```

We often speak of double-precision floating-point numbers as having about 16 decimal digits. The 52-bit significand is paired with a sign bit and 11 binary bits to represent the exponent $n$ in {eq}`floatpoint`, for a total of 64 binary bits per floating-point number.

(demo-float-double)=
::::{prf:example} Floating-point representation
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{include} julia/double.ipynb
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-float-double-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-float-double-python
:::
```` 
`````
::::

Our theoretical description of $\float$ did not place limits on the exponent, but in double precision its range is limited to $-1022\le n \le 1023$. Thus, the largest number is just short of $2^{1024}\approx 2\times 10^{308}$, which is enough in most applications. Results that should be larger are said to *overflow* and will actually result in the value `Inf`. Similarly, the smallest positive number is $2^{-1022}\approx 2\times 10^{-308}$, and smaller values are said to *underflow* to zero.[^denormalized]

[^denormalized]: Actually, there are some still-smaller *denormalized* numbers that have less precision, but we won't use that level of detail.

Note the crucial difference between $\macheps=2^{-52}$, which is the distance between 1 and the next larger double-precision number, and $2^{-1022}$, which is the smallest positive double-precision number. The former has to do with relative precision, while the latter is about absolute precision. Getting close to zero always requires a shift in thinking to absolute precision because any finite error is infinite relative to zero.

```{index} NaN
```

One more double-precision value is worthy of note: **NaN**, which stands for **Not a Number**. It is the result of an undefined arithmetic operation such as 0/0.

## Floating-point arithmetic

```{index} floating-point numbers
```

Computer arithmetic is performed on floating-point numbers and returns floating-point results. We assume the existence of machine-analog operations for real functions such as $+$, $-$, $\times$, $/$, $\sqrt{\quad}$, and so on. Without getting into the details, we will suppose that each elementary machine operation creates a floating-point result whose relative error is bounded by $\macheps$. For example, if $x$ and $y$ are in $\float$, then for machine addition $\oplus$ we have the bound

```{math}
\frac{ |(x \oplus y)-(x+y)| }{ |x+y| } \le \macheps.
```
Hence the relative error in arithmetic is essentially the same as for the floating-point representation itself. However, playing by these rules can lead to disturbing results.

(demo-float-arithmetic)=
::::{prf:example} Floating-point arithmetic oddity
`````{tab-set} 
````{tab-item} Julia
:sync: julia
:::{include} julia/arithmetic.ipynb
:::
```` 

````{tab-item} MATLAB
:sync: matlab
:::{embed} #demo-float-arithmetic-matlab
:::
```` 

````{tab-item} Python
:sync: python
:::{embed} #demo-float-arithmetic-python
:::
```` 
`````
::::

There are two ways to look at {numref}`Demo %s <demo-float-arithmetic>`. On one hand, its two versions of the result differ by less than $1.2\times 10^{-16}$, which is very small — not just in everyday terms, but with respect to the operands, which are all close to 1 in absolute value. On the other hand, the difference is as large as the exact result itself! We formalize and generalize this observation in the next section. In the meantime, keep in mind that exactness cannot be taken for granted in floating-point computation. 

::::{prf:observation}
We should not expect that two mathematically equivalent results will be equal when computed in floating point, only that they be relatively close together.
::::

## Exercises

```{note}
Exercises marked with ✍ are intended to be done by hand or with the aid of a simple calculator. Exercises marked with ⌨ are intended to be solved using a computer.
```
(problem-float-d4)=
1. ✍ Consider a floating-point set $\float$ defined by {eq}`floatpoint` and {eq}`mantissa` with $d=4$.
  
    **(a)** How many elements of $\float$ are there in the real interval $[1/2,4]$, including the endpoints?

    **(b)** What is the element of $\float$ closest to the real number $1/10$? (Hint: Find the interval $[2^n,2^{n+1})$ that contains $1/10$, then enumerate all the candidates in $\float$.)

    **(c)** What is the smallest positive integer not in $\float$? (Hint: For what value of the exponent does the spacing between floating-point numbers become larger than 1?)

(problem-fp-fprelative)=
2. ✍ Prove that {eq}`fpbound` is equivalent to {eq}`fpboundalt`. This means showing first that {eq}`fpbound` implies {eq}`fpboundalt`, and then separately that {eq}`fpboundalt` implies {eq}`fpbound`.

3. ⌨ There are much better rational approximations to $\pi$ than 22/7 as used in {numref}`Demo {number} <demo-float-accuracy>`. For each one below, find its absolute and relative accuracy, and (rounding down to an integer) the number of accurate digits. 

    **(a)** 355/113

    **(b)** 103638/32989

    ```{index} IEEE 754
    ```

4. ✍ IEEE 754 **single precision** specifies that 23 binary bits are used for the value $f$ in the significand $1+f$ in {eq}`mantissa`. Because they need less storage space and can be operated on more quickly than double-precision values, single-precision values can be useful in low-precision applications. (They are supported as type `Float32` in Julia.)

    **(a)** In base-10 terms, what is the first single-precision number greater than $1$ in this system?

    **(b)** What is the smallest positive integer that is not a single-precision number? (See the hint to Exercise 1.)

5. ⌨ Julia defines a function `nextfloat` that gives the next-larger floating-point value of a given number. What is the next float past `floatmax()`? What is the next float past `-Inf`? 
