(chapter-leastsq)=

# 3. Overdetermined linear systems

```{index} Han Solo, The Empire Strikes Back
```

```{epigraph}
I must have hit pretty close to the mark to get her all riled up like that, huh, kid? 

— Han Solo, *The Empire Strikes Back*
```

So far we have considered $\mathbf{A}\mathbf{x}=\mathbf{b}$ only when $\mathbf{A}$ is a square matrix. In this chapter we consider how to interpret and solve the problem for an $m\times n$ matrix where $m>n$—and in practice, $m$ is often *much* larger than $n$. This is called an *overdetermined* linear system because, in general, the system has more equations to satisfy than the variables allow. The complementary *underdetermined* case $m<n$ turns up less frequently and will not be considered in this book.

Since we cannot solve all the system's equations, we need to define what the "best possible" answer is. There are multiple useful options, but the most important version of the overdetermined problem occurs using the *least squares*—the sum of the squares of the equation residuals is minimized. This is far from an arbitrary choice. Mathematically, we recognize the sum-of-squares as a vector 2-norm and therefore tied to inner products; physically, the 2-norm may coincide with energy, which is often minimized by natural systems; and statistically, least squares leads to the estimates of maximum likelihood for certain models. Furthermore, the solution of the least-squares problem requires only linear algebra and is about as easily to compute as in the square case.

The linear least-squares problem serves as our introduction to the vast field of *optimization*. It is one of the simplest problems of this type. We will see an extension to a nonlinear version in the next chapter.
