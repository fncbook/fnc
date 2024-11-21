---
jupytext:
  cell_metadata_filter: -all
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Julia 1.7.1
  language: julia
  name: julia-fast
---
```{code-cell}
:tags: [remove-cell]
using FundamentalsNumericalComputation
FNC.init_format()
```

(section-nonlineqn-quasinewton)=
# Quasi-Newton methods

Newton's method is a foundation for algorithms to solve equations and minimize quantities. But it is not ideal in its straightforward or pure form. Specifically, its least appealing features are the programming nuisance and computational expense of evaluating the Jacobian matrix, and the tendency of the iteration to diverge from many starting points. There are different **quasi-Newton methods** that modify the basic idea in an attempt to overcome these issues.

## Jacobian by finite differences

```{index} Jacobian matrix, finite differences
```

In the scalar case, we found an easy alternative to a direct evaluation of the derivative. In retrospect, we may interpret the secant formula {eq}`secant` as the Newton formula {eq}`newton` with $f'(x_k)$ replaced by the difference quotient

```{math}
  :label: secantfd
  \frac{f(x_k)-f(x_{k-1})}{x_k-x_{k-1}}.
```

If the sequence of $x_k$ values converges to a root $r$, then this quotient converges to $f'(r)$.

In the system case, replacing the Jacobian evaluation is more complicated: derivatives are needed with respect to $n$ variables, not just one. From {eq}`jacobian`, we note that the $j$th column of the Jacobian is

```{math}
  \mathbf{J}(\mathbf{x}) \mathbf{e}_j =
  \begin{bmatrix}
    \frac{\partial{f_1}}{\partial x_j} \\[2mm] \frac{\partial{f_2}}{\partial x_j}
    \\ \vdots \\ \frac{\partial{f_n}}{\partial x_j}
  \end{bmatrix}.
```

(As always, $\mathbf{e}_j$ represents the $j$th column of the identity matrix, here in $n$ dimensions.) Inspired by {eq}`secantfd`, we can replace the differentiation with a quotient involving a change in only $x_j$ while the other variables remain fixed:

```{math}
  :label: jacobianfd
  \mathbf{J}(\mathbf{x}) \mathbf{e}_j \approx
  \frac{\mathbf{f}(\mathbf{x}+\delta \mathbf{e}_j) - \mathbf{f}(\mathbf{x})}{\delta}, \qquad j=1,\ldots,n.
```

```{index} machine epsilon
```

For reasons explained in Chapter 5, $\delta$ is usually chosen close to $\sqrt{\epsilon}$, where $\epsilon$ represents the expected noise or uncertainty level in evaluation of $\mathbf{f}$. If the only source of noise is floating-point roundoff, then $\delta \approx \sqrt{\epsilon_\text{mach}}$.

The finite-difference formula {eq}`jacobianfd` is implemented by {numref}`Function {number} <function-fdjac>`.

(function-fdjac)=
````{prf:function} fdjac
**Finite-difference approximation of a Jacobian**

```{code-block} julia
:lineno-start: 1
"""
    fdjac(f,x₀[,y₀])

Compute a finite-difference approximation of the Jacobian matrix for
`f` at `x₀`, where `y₀`=`f(x₀)` may be given.
"""
function fdjac(f,x₀,y₀=f(x₀))
    δ = sqrt(eps())*max(norm(x₀),1)   # FD step size
    m,n = length(y₀),length(x₀)
    if n==1
        J = (f(x₀+δ) - y₀) / δ
    else
        J = zeros(m,n)
        x = copy(x₀)
        for j in 1:n
            x[j] += δ
            J[:,j] = (f(x) - y₀) / δ
            x[j] -= δ
        end
    end
    return J
end
```
````

::::{admonition} About the code
:class: dropdown
{numref}`Function {number} <function-fdjac>` is written to accept the case where $\mathbf{f}$ maps $n$ variables to $m$ values with $m\neq n$, in anticipation of {numref}`section-nonlineqn-nlsq`.

Note that a default value is given for the third argument `y₀`, and it refers to earlier arguments in the list. The reason is that in some contexts, the caller of `fdjac` may have already computed `y₀` and can supply it without computational cost, while in other contexts, it must be computed fresh. The configuration here adapts to either situation.
::::

## Broyden's update

The finite-difference Jacobian is easy to conceive and use. But, as you can see from {eq}`jacobianfd`, it requires $n$ additional evaluations of the system function at each iteration, which can be unacceptably slow in some applications. Conceptually these function evaluations seem especially wasteful given that the root estimates, and thus presumably the Jacobian matrix, are supposed to change little as the iteration converges. This is a good time to step in with the principle of approximate approximation, which suggests looking for a shortcut in the form of a cheap-but-good-enough way to update the Jacobian from one iteration to the next.

Recall that the Newton iteration is derived by solving the linear model implied by {eq}`multitaylor`:

```{math}
  \mathbf{f}(\mathbf{x}_{k+1}) \approx \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)\,(\mathbf{x}_{k+1}-\mathbf{x}_k) = \boldsymbol{0}.
```

Let $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$  be the Newton step. Let $\mathbf{y}_k=\mathbf{f}(\mathbf{x}_k)$, and now we replace $\mathbf{J}(\mathbf{x}_k)$ by a matrix $\mathbf{A}_{k}$ that is meant to approximate the Jacobian. Hence the Newton step is considered to be defined, as in {numref}`Algorithm {number} <algorithm-nonlineqn-newtonsys>`, by

```{math}
  :label: quasinewton-step
  \mathbf{A}_k \mathbf{s}_k = -\mathbf{y}_k.
```

Once $\mathbf{x}_{k+1}$ is obtained, we should update the approximate Jacobian to a new $\mathbf{A}_{k+1}$. If we think one-dimensionally for a moment, the secant method would assume that $A_{k+1}=(f_{k+1}-f_k)/(x_{k+1}-x_k)$. It's not easy to generalize a fraction to vectors, but we can do it if we instead write it as

```{math}
  \mathbf{y}_{k+1}-\mathbf{y}_k = \mathbf{A}_{k+1} (\mathbf{x}_{k+1}-\mathbf{x}_k) = \mathbf{A}_{k+1} \mathbf{s}_k.
```

This is used to justify the following requirement:

```{math}
  :label: secantsys
  \mathbf{A}_{k+1} \mathbf{s}_k = \mathbf{y}_{k+1}-\mathbf{y}_k.
```

```{index} ! Broyden update
```

This isn't enough to uniquely determine $\mathbf{A}_{k+1}$. However, if we also require that $\mathbf{A}_{k+1}-\mathbf{A}_k$ is a matrix of rank 1, then one arrives at the following.

::::{prf:definition} Broyden update formula
Using the definitions above,

```{math}
  :label: broyden
  \mathbf{A}_{k+1} = \mathbf{A}_k + \frac{1}{\mathbf{s}_k^T \mathbf{s}_k}(\mathbf{y}_{k+1} - \mathbf{y}_k -\mathbf{A}_k \mathbf{s}_k)\, \mathbf{s}_k^T.
```
::::

```{index} outer product
```

Observe that $\mathbf{A}_{k+1}-\mathbf{A}_k$ is proportional to the outer product of two vectors, and that computing it requires no extra evaluations of $\mathbf{f}$. Remarkably, under reasonable assumptions, the sequence of $\mathbf{x}_k$ resulting when Broyden updates are used converges superlinearly, even though the matrices $\mathbf{A}_k$ do not necessarily converge to the Jacobian of $\mathbf{f}$. 

In practice, one typically uses finite differences to initialize the Jacobian at iteration $k=1$. If for some $k$ the step computed by the update formula fails to make enough improvement in the residual, then $\mathbf{A}_k$ is reinitialized by finite differences and the step is recalculated.

## Levenberg's method

The most difficult part of many rootfinding problems is finding a starting point that will lead to convergence. The linear model implicitly constructed during a Newton iteration—whether we use an exact, finite-difference, or iteratively updated Jacobian matrix—becomes increasingly inaccurate as one ventures farther from the most recent root estimate, eventually failing to resemble the exact function much at all. 

Although one could imagine trying to do a detailed accuracy analysis of each linear model as we go, in practice simple strategies are valuable here. Suppose, after computing the step suggested by the linear model, we ask a binary question: Would taking that step improve our situation? Since we are trying to find a root of $\mathbf{f}$, we have a quantitative way to pose this question: Does the backward error $\|\mathbf{f}\|$ decrease? If not, we should reject the step and find an alternative.

There are several ways to find alternatives to the standard step, but we will consider just one of them, based on the parameterized equation 

```{math}
  :label: levenberg
  (\mathbf{A}_k^T \mathbf{A}_k + \lambda \mathbf{I})\,\mathbf{s}_k = -\mathbf{A}_k^T \mathbf{f}_k.
```

```{index} ! Levenberg's method
```

(algorithm-nonlineqn-levenberg)=
::::{prf:algorithm} Levenberg's method
Given $\mathbf{f}$, a starting value $\mathbf{x}_1$, and a scalar $\lambda$, for each $k=1,2,3,\ldots$

1. Compute $\mathbf{y}_k = \mathbf{f}(\mathbf{x}_k)$, and let $\mathbf{A}_k$ be an exact or approximate Jacobian matrix.
2. Solve the linear system {eq}`levenberg` for $\mathbf{s}_k$.
3. Let $\hat{\mathbf{x}} = \mathbf{x}_k + \mathbf{s}_k$.
4. If the residual is reduced at $\hat{\mathbf{x}}$, then let $\mathbf{x}_{k+1}=\hat{\mathbf{x}}$.
5. Update $\lambda$ and update $\mathbf{A}_k$ to $\mathbf{A}_{k+1}$.
::::

Some justification of {eq}`levenberg` comes from considering extreme cases for $\lambda$. If $\lambda=0$, then

```{math}
  \mathbf{A}_k^T \mathbf{A}_k \mathbf{s}_k = -\mathbf{A}_k^T \mathbf{f}_k,
```

which is equivalent to the definition of the usual linear model (i.e., Newton or quasi-Newton) step {eq}`quasinewton-step`. On the other hand, as $\lambda\to\infty$, Equation {eq}`levenberg` approaches

```{math}
  :label: steepest
  \lambda \mathbf{s}_k = - \mathbf{A}_k^T \mathbf{f}_k.
```

To interpret this equation, define the scalar residual function 

$$
\phi(\mathbf{x})=\mathbf{f}(\mathbf{x})^T\mathbf{f}(\mathbf{x}) = \|\mathbf{f}(\mathbf{x})\|^2.
$$

Finding a root of $\mathbf{f}$ is equivalent to minimizing $\phi$. A calculation shows that the gradient of $\phi$ is

```{math}
   :label: nlsgradient
   \nabla \phi(\mathbf{x}) = 2 \mathbf{J}(\mathbf{x})^T \mathbf{f}(\mathbf{x}).
```

```{index} steepest descent
```

Hence if $\mathbf{A}_k=\mathbf{J}(\mathbf{x}_k)$, then $\mathbf{s}_k$ from {eq}`steepest` is in the opposite direction from the gradient vector. In vector calculus you learn that this direction is the one of most rapid decrease or **steepest descent**. A small enough step in this direction is guaranteed in all but pathological cases to decrease $\phi$, which is exactly what we want from a backup plan.

In effect, the $\lambda$ parameter in {eq}`levenberg` allows a smooth transition between the pure Newton step, for which convergence is very rapid near a root, and a small step in the gradient descent direction, which guarantees progress for the iteration when we are far from a root. 

## Implementation

To a large extent, the incorporation of finite differences, Jacobian updates, and Levenberg step are independent decisions. {numref}`Function {number} <function-levenberg>` shows how they might be combined. This function is one of the most logically complex we have encountered so far.

Each pass through the loop starts by using {eq}`levenberg` to propose a step $\mathbf{s}_k$. The function then asks whether using this step would decrease the value of $\|\mathbf{f}\|$ from its present value. If so, we accept the new root estimate, we decrease $\lambda$ in order to get more Newton-like (since things have gone well), and we apply the Broyden formula to get a cheap update of the Jacobian. If the proposed step is not successful, we increase $\lambda$ to get more gradient-like (since we just failed) and, if the current Jacobian was the result of a cheap update, use finite differences to reevaluate it.  

(function-levenberg)=
````{prf:function} levenberg
**Quasi-Newton method for nonlinear systems**

```{code-block} julia
:lineno-start: 1
"""
    levenberg(f,x₁[;maxiter,ftol,xtol])

Use Levenberg's quasi-Newton iteration to find a root of the system
`f` starting from `x₁`. Returns the history of root estimates 
as a vector of vectors.

The optional keyword parameters set the maximum number of iterations
and the stopping tolerance for values of `f` and changes in `x`.

"""
function levenberg(f,x₁;maxiter=40,ftol=1e-12,xtol=1e-12)
    x = [float(x₁)]
    yₖ = f(x₁)
    k = 1;  s = Inf;
    A = fdjac(f,x[k],yₖ)   # start with FD Jacobian
    jac_is_new = true

    λ = 10;
    while (norm(s) > xtol) && (norm(yₖ) > ftol)
        # Compute the proposed step.
        B = A'*A + λ*I
        z = A'*yₖ
        s = -(B\z)
        
        x̂ = x[k] + s
        ŷ = f(x̂)

        # Do we accept the result?
        if norm(ŷ) < norm(yₖ)    # accept
            λ = λ/10   # get closer to Newton
            # Broyden update of the Jacobian.
            A += (ŷ-yₖ-A*s)*(s'/(s'*s))
            jac_is_new = false
            
            push!(x,x̂)
            yₖ = ŷ
            k += 1
        else                       # don't accept
            # Get closer to gradient descent.
            λ = 4λ
            # Re-initialize the Jacobian if it's out of date.
            if !jac_is_new
                A = fdjac(f,x[k],yₖ)
                jac_is_new = true
            end
        end

        if k==maxiter
            @warn "Maximum number of iterations reached."
            break
        end
        
    end
    return x
end
```
````

In some cases our simple logic in {numref}`Function {number} <function-levenberg>` can make $\lambda$ oscillate between small and large values; several better but more complicated strategies for controlling $\lambda$ are known. In addition, the linear system {eq}`levenberg` is usually modified to get the well-known **Levenberg–Marquardt** algorithm, which does a superior job in some problems as $\lambda\to \infty$.

(demo-quasi-levenberg)=
```{prf:example}
```





To solve a nonlinear system, we need to code only the function defining the system, and not its Jacobian.

```{code-cell}
function f(x)  
    [
      exp(x[2]-x[1]) - 2,
      x[1]*x[2] + x[3],
      x[2]*x[3] + x[1]^2 - x[2]
    ]
end;
```

In all other respects usage is the same as for the `newtonsys` function.

```{code-cell}
x₁ = [0.,0.,0.]   
x = FNC.levenberg(f,x₁)
```

It's always a good idea to check the accuracy of the root, by measuring the residual (backward error).

```{code-cell}
r = x[end]
println("backward error = $(norm(f(r)))")
```

Looking at the convergence in norm, we find a convergence rate between linear and quadratic, like with the secant method.

```{code-cell}
logerr = [ log( norm(x[k]-r) ) for k in 1:length(x)-1 ]
[ logerr[k+1]/logerr[k] for k in 1:length(logerr)-1 ]
```




## Exercises

1. ⌨ (Variation on [Exercise 4.5.2](problem-newtonsys-spherepotential).) Two curves in the $(u,v)$ plane are defined implicitly by the equations $u\log u + v \log v = -0.3$ and $u^4 + v^2 = 1$.
  
    **(a)** ✍ Write the intersection of these curves in the form $\mathbf{f}(\mathbf{x}) = \boldsymbol{0}$ for two-dimensional $\mathbf{f}$ and $\mathbf{x}$.

    **(b)** ⌨ Use {numref}`Function {number} <function-levenberg>` to find an intersection point starting from $u=1$, $v=0.1$.

    **(d)** ⌨ Use {numref}`Function {number} <function-levenberg>` to find an intersection point starting from $u=0.1$, $v=1$.

2. ⌨ (Variation on [Exercise 4.5.4](problem-newtonsys-orbitintersect)) Two elliptical orbits $(x_1(s),y_1(s))$ and $(x_2(t),y_2(t))$ are described by the equations
  
    ```{math}
    \begin{bmatrix}
      x_1(t) \\ y_1(t)
    \end{bmatrix}
    =
    \begin{bmatrix}
      -5+10\cos(t) \\ 6\sin(t)
    \end{bmatrix}, \qquad
    \begin{bmatrix}
      x_2(t)\\y_2(t)
    \end{bmatrix} =
    \begin{bmatrix}
      8\cos(t) \\ 1+12\sin(t)
    \end{bmatrix},
    ```

    where $t$ represents time.

    **(a)** ✍ Write out a $2\times 2$ nonlinear system of equations that describes an intersection of these orbits. (Note: An intersection is not the same as a collision—they don't have to occupy the same point at the same time.)

    **(b)** ⌨ Use {numref}`Function {number} <function-levenberg>` to find all of the unique intersections.

3. ⌨  (Variation on [Exercise 4.5.5](problem-newtonsys-ellipsemin)) Suppose one wants to find the points on the ellipsoid $x^2/25 + y^2/16 + z^2/9 = 1$ that are closest to and farthest from the point $(5,4,3)$. The method of Lagrange multipliers implies that any such point satisfies
  
    ```{math}
    \begin{split}
        x-5 &= \frac{\lambda x}{25}, \\[1mm]
        y-4 &= \frac{\lambda y}{16}, \\[1mm]
        z-3 &= \frac{\lambda z}{9}, \\[1mm]
        1 &=  \frac{1}{25}x^2 + \frac{1}{16}y^2 + \frac{1}{9}z^2
    \end{split}
    ```

    for an unknown value of $\lambda$.
  
    **(a)** Write out this system in the form $\mathbf{f}(\mathbf{u}) = \boldsymbol{0}$. (Note that the system has four variables to go with the four equations.)

    **(b)** Use {numref}`Function {number} <function-levenberg>` with different initial guesses to find the two roots of this system. Which is the closest point to $(5,4,3)$, and which is the farthest?

 
4. ✍ The Broyden update formula {eq}`broyden` is just one instance of so-called rank-1 updating. Verify the  *Sherman–Morrison formula*,
  
    ```{math}
    (\mathbf{A}+\mathbf{u}\mathbf{v}^T)^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\frac{\mathbf{u}\mathbf{v}^T}{1+\mathbf{v}^T\mathbf{A}^{-1}\mathbf{u}}\mathbf{A}^{-1},
    ```

    which is valid whenever $\mathbf{A}$ is invertible and the denominator above is nonzero. (Hint: Show that $\mathbf{A}+\mathbf{u}\mathbf{v}^T$ times the matrix above simplifies to the identity matrix.)

5. ✍ Derive Equation {eq}`nlsgradient`.
  
6. ⌨ (See also [Exercise 4.5.11](problem-newtonsys-byhand).) Suppose that

    ```{math}
    \mathbf{f}(\mathbf{x}) =
    \begin{bmatrix}
      x_1x_2+x_2^2-1 \\[1mm] x_1x_2^3 + x_1^2x_2^2 + 1
    \end{bmatrix}.
    ```

    Let $\mathbf{x}_1=[-2,1]^T$ and let $\mathbf{A}_1=\mathbf{J}(\mathbf{x}_1)$ be the exact Jacobian.
  
    **(a)** Solve {eq}`levenberg` for $\mathbf{s}_1$ with $\lambda=0$; this is the "pure" Newton step. Show numerically that $\|\mathbf{f}(\mathbf{x}_1+\mathbf{s}_1)\| > \|\mathbf{f}(\mathbf{x}_1)\|$. (Thus, the Newton step made us go to a point seemingly farther from a root than where we started.)

    **(b)** Now repeat part (a) with $\lambda=0.01j$ for $j=1,2,3,\ldots.$ What is the smallest value of $j$ such that $\|\mathbf{f}(\mathbf{x}_1+\mathbf{s}_1)\| < \|\mathbf{f}(\mathbf{x}_1)\|$?

7. ✍ Show that Equation {eq}`levenberg` is equivalent to the linear least-squares problem
  
    ```{math}
    \min_{\mathbf{v}} \Bigl(  \bigl\|\mathbf{A}_k\mathbf{v} + \mathbf{f}_k\bigr\|_2^2 +
    \lambda^2 \bigl\| \mathbf{v} \bigr\|_2^2 \Bigr).
    ```

    (Hint: Express the minimized quantity using block matrix notation, such that {eq}`levenberg` becomes the normal equations for it.)

    Thus, another interpretation of Levenberg's method is that it is the Newton step plus a penalty, weighted by $\lambda$, for taking large steps.


