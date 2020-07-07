# Quasi-Newton methods

Newton's method is a foundation for algorithms to solve equations and minimize quantities. But it is not ideal in its plain or "pure" form. Instead there are different {term}`quasi-Newton methods` that attempt to overcome two serious issues: the programming nuisance and computational expense of evaluating the Jacobian matrix, and the tendency of the iteration to diverge for many starting points.

## Jacobian by finite differences

In the scalar case, we found an easy alternative to a direct evaluation of the derivative. Specifically, we may interpret the secant formula {eq}`secant` as the Newton formula {eq}`newton` with $f'(x_k)$ replaced by the quotient

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

For reasons explained in the next chapter, $\delta$ is usually chosen close to $\sqrt{\epsilon}$, where $\epsilon$ represents the expected noise level in evaluation of $\mathbf{f}$. If the only source of noise is floating-point roundoff, then $\delta=\sqrt{\epsilon_\text{mach}}$.

The finite-difference formula {eq}`jacobianfd` is implemented by the short code {ref}`function-fdjac`. (The code is written to accept the case where $\mathbf{f}$ maps $n$ variables to $m$ values with $m\neq n$, in anticipation of \secref{nl-least-sq}.)

(function-fdjac)=

````{proof:function} fdjac
**Finite-difference approximation of a Jacobian.**

```{code-block} julia
:lineno-start: 1
"""
fdjac(f,x0,y0)

Compute a finite-difference approximation of the Jacobian matrix for
`f` at `x0`, where `y0`=`f(x0)` is given.
"""
function fdjac(f,x0,y0)

delta = sqrt(eps())   # FD step size
m,n = length(y0),length(x0)
if n==1
    J = (f(x0+delta) - y0) / delta
else
    J = zeros(m,n)
    In = I(n)
    for j = 1:n
        J[:,j] = (f(x0 + delta*In[:,j]) - y0) / delta
    end
end

return J
end
```
````

## Broyden's update

The finite-difference Jacobian is easy to conceive and use. But, as you can see from {eq}`jacobianfd`, it requires $n$ additional evaluations of the system function at each iteration, which can be unacceptably slow in some applications. Conceptually these function evaluations seem especially wasteful given that the root estimates, and thus presumably the Jacobian matrix, are supposed to change little as the iteration converges. This is a good time to step in with the principle of approximate approximation, which suggests looking for a shortcut in the form of a cheap-but-good-enough way to update the Jacobian from one iteration to the next.

Recall that the Newton iteration is derived by solving the linear model implied by {eq}`multitaylor`:

```{math}
  \boldsymbol{0} \approx \mathbf{f}(\mathbf{x}_{k+1}) \approx \mathbf{f}(\mathbf{x}_k) + \mathbf{J}(\mathbf{x}_k)\,(\mathbf{x}_{k+1}-\mathbf{x}_k).
```

Let $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$  be the Newton step. We will make the notation simpler via $\mathbf{f}_k=\mathbf{f}(\mathbf{x}_k)$, and now we replace $\mathbf{J}(\mathbf{x}_k)$ by a matrix $\mathbf{A}_{k}$ that is meant to approximate the Jacobian. Hence the Newton step is considered to be defined by

```{math}
  :label: quasinewton-step
  \mathbf{A}_k \mathbf{s}_k = -\mathbf{f}_k.
```

This equation gets us from $\mathbf{x}_k$ to $\mathbf{x}_{k+1}$. To continue the iteration, we want to update the approximate Jacobian to $\mathbf{A}_{k+1}$. If we think one-dimensionally for a moment, the secant method would assume that $A_{k+1}=(f_{k+1}-f_k)/(x_{k+1}-x_k)$. It's not easy to generalize a fraction to vectors, but we can do it if we instead write it as

```{math}
  \mathbf{f}_{k+1}-\mathbf{f}_k = \mathbf{A}_{k+1} (\mathbf{x}_{k+1}-\mathbf{x}_k) = \mathbf{A}_{k+1} \mathbf{s}_k.
```

This is used to justify the following requirement:

```{math}
  :label: secantsys
  \mathbf{A}_{k+1} \mathbf{s}_k = \mathbf{f}_{k+1}-\mathbf{f}_k.
```

This isn't enough to uniquely determine $\mathbf{A}_{k+1}$. However, if we also require that $\mathbf{A}_{k+1}-\mathbf{A}_k$ is a matrix of rank one, then one arrives at the **Broyden update formula** \index{Broyden update} \index{quasi-Newton methods!Broyden update}

```{math}
  :label: broyden
  \mathbf{A}_{k+1} = \mathbf{A}_k + \frac{1}{\mathbf{s}_k^T \mathbf{s}_k}(\mathbf{f}_{k+1} - \mathbf{f}_k -\mathbf{A}_k \mathbf{s}_k)\, \mathbf{s}_k^T.
```

Observe that $\mathbf{A}_{k+1}-\mathbf{A}_k$, being proportional to the outer product of two vectors, is indeed a rank-one matrix, and that computing it requires no extra evaluations of $\mathbf{f}$. Remarkably, under reasonable assumptions the sequence of $\mathbf{x}_k$ so defined converges superlinearly, even though the matrices $\mathbf{A}_k$ do not necessarily converge to the Jacobian of $\mathbf{f}$. In practice one typically uses finite differences to initialize the Jacobian at iteration $k=1$. If the step computed by the update formula improves the solution, it is accepted and the iteration continues, but if the update formula fails to give a good result, the matrix $\mathbf{A}_k$ is reinitialized by finite differences and the step is recalculated.

## Levenberg's method

The most difficult part of many rootfinding problems is finding a starting point that will lead to convergence. The linear model implicitly constructed during a Newton iteration---whether we use an exact, finite-difference, or iteratively updated Jacobian matrix---becomes increasingly inaccurate as one ventures farther from the most recent root estimate, eventually failing to resemble the exact function much at all. Although one could imagine trying to do a detailed accuracy analysis of each linear model as we go, in practice simple strategies are valuable here. Suppose, after computing the step suggested by the linear model, we ask a binary question: Would taking that step improve our situation? Since we are trying to find a root of $\mathbf{f}$, we have a quantitative way to pose this question: Does the backward error $\|\mathbf{f}\|$ decrease? If not, we should reject the step and find an alternative.

There are several ways to find alternatives to the standard step, but we will consider just one of them. Let $\mathbf{A}_k$ be the (exact or approximate) Jacobian matrix for iteration number $k$. **Levenberg's method** introduces a positive parameter $\lambda$ into the calculation of the next step: define

```{math}
  :label: levenberg
  (\mathbf{A}_k^T \mathbf{A}_k + \lambda \mathbf{I})\,\mathbf{s}_k = -\mathbf{A}_k^T \mathbf{f}_k,
```

where $\mathbf{x}_{k+1}=\mathbf{x}_k+\mathbf{s}_k$. Some justification of {eq}`levenberg` comes from considering extreme cases for $\lambda$. If $\lambda=0$, then

```{math}
  \mathbf{A}_k^T \mathbf{A}_k \mathbf{s}_k = -\mathbf{A}_k^T \mathbf{f}_k,
```

which is equivalent to the definition of the usual linear model (i.e., Newton or quasi-Newton) step {eq}`quasinewton-step`. On the other hand, as $\lambda\to\infty$, equation {eq}`levenberg` is increasingly close to

```{math}
  :label: steepest
  \lambda \mathbf{s}_k = - \mathbf{A}_k^T \mathbf{f}_k.
```

To interpret this equation, define the scalar residual function $\phi(\mathbf{x})=\mathbf{f}(\mathbf{x})^T\mathbf{f}(\mathbf{x}) = \|\mathbf{f}(\mathbf{x})\|^2$. Finding a root of $\mathbf{f}$ is equivalent to minimizing $\phi$. A calculation shows that the gradient of $\phi$ is

```{math}
   :label: nlsgradient
   \nabla \phi(\mathbf{x}) = 2 \mathbf{J}(\mathbf{x})^T \mathbf{f}(\mathbf{x}).
```

Hence if $\mathbf{A}_k=\mathbf{J}(\mathbf{x}_k)$, then $\mathbf{s}_k$ from {eq}`steepest` is in the opposite direction from the gradient vector. In vector calculus you learn that this direction is the one of most rapid decrease; for this reason {eq}`steepest` is called the **steepest descent direction.** \index{steepest descent} A small enough step in this direction is guaranteed (in all but pathological cases) to decrease $\phi$, which is exactly what we want from a backup plan.

In summary, the $\lambda$ parameter in {eq}`levenberg` allows a smooth transition between the pure Newton step, for which convergence is very rapid (i.e., near a root), and a small step in the descent direction, which guarantees some progress for the iteration when we are far from a root. To make the Levenberg step computation into an algorithm, we will combine it with an accept/reject strategy as described above.

## Implementation

````{sidebar} Demo
:class: demo
{doc}`demos/quasi-levenberg`
````

To a large extent the incorporation of finite differences, Jacobian updates, and Levenberg step are independent. {ref}`function-levenberg` shows how they might be combined. This function is one of the most logically complex we have encountered so far.

First observe that we have introduced a MATLAB convenience feature: *optional input parameters*. The keyword "nargin" inside a function evaluates to the number of input arguments that were provided by the caller. In this case we supply a default value for the third argument, "tol", if none was specified by the caller. Most modern computing languages have  analogous mechanisms for accepting a variable number of input parameters. A common use case is what we have done here, allowing the optional override of a default setting.

Each pass through the loop starts by using {eq}`levenberg` to propose a step $\mathbf{s}_k$. The algorithm then asks whether using this step would decrease the value of $\|\mathbf{f}\|$ from its present value. If so, $\mathbf{x}_k+\mathbf{s}_k$ is the new root estimate; since the iteration is going well, we decrease $\lambda$ (i.e., get more Newton-like) and apply the Broyden formula to get a fast update of the Jacobian. If the proposed step is not successful, we increase $\lambda$ (more descent-like) and, if the current Jacobian was the result of a cheap update, use finite differences to reevaluate it. Whether or not the step was accepted, everything is set up for the next loop iteration. Note that the loop iterations no longer correspond to the (quasi-)Newton iteration counter $k$.

Finally, we draw attention to lines 50--52. Rather than issuing a warning if the number of iterations is too large, we do it when the final residual is fairly large. This is done to avoid silently returning a value that is easily seen to be incorrect.

(function-levenberg)=

````{proof:function} levenberg
**Quasi-Newton method for nonlinear systems.**

```{code-block} julia
:lineno-start: 1
"""
levenberg(f,x1,tol)

Use Levenberg's quasi-Newton iteration to find a root of the system
`f`, starting from `x1`, with `tol` as the stopping tolerance in
both step size and residual norm. Returns root estimates as a
matrix, one estimate per column.
"""
function levenberg(f,x1,tol=1e-12)

# Operating parameters.
ftol = tol;  xtol = tol;  maxiter = 40;

x = zeros(length(x1),maxiter)
x = [float(x1)]
fk = f(x1)
k = 1;  s = Inf;
Ak = fdjac(f,x1,fk)   # start with FD Jacobian
jac_is_new = true

lambda = 10;
while (norm(s) > xtol) && (norm(fk) > ftol) && (k < maxiter)
    # Compute the proposed step.
    B = Ak'*Ak + lambda*I
    z = Ak'*fk
    s = -(B\z)

    xnew = x[k] + s
    fnew = f(xnew)

    # Do we accept the result?
    if norm(fnew) < norm(fk)    # accept
        y = fnew - fk
        push!(x,xnew)
        fk = fnew
        k += 1

        lambda = lambda/10   # get closer to Newton
        # Broyden update of the Jacobian.
        Ak = Ak + (y-Ak*s)*(s'/(s'*s))
        jac_is_new = false
    else                       # don't accept
        # Get closer to steepest descent.
        lambda = 4lambda
        # Re-initialize the Jacobian if it's out of date.
        if !jac_is_new
            Ak = fdjac(f,x[k],fk)
            jac_is_new = true
        end
    end
end

if (norm(fk) > 1e-3)
    @warn "Iteration did not find a root."
end

return x
end
```
````

<!-- ````{sidebar} Demo
:class: demo
{doc}`demos/quasi-MM`
```` 
-->

In some cases our simple logic in {ref}`function-levenberg` can make $\lambda$ oscillate between small and large values; several better but more complicated strategies for controlling $\lambda$ are known. In addition, the linear system {eq}`levenberg` is usually modified to get the well-known **Levenberg--Marquardt** algorithm, which does a superior job in some problems as $\lambda\to \infty$.

## Exercises

1. ⌨ (variation of [earlier problem](problem-spherepotential)) Two curves in the $(u,v)$ plane are defined implicitly by the equations $u\log u + v \log v = -0.3$ and $u^4 + v^2 = 1$.
  
    **(a)** ✍ Write the intersection of these curves in the form $\mathbf{f}(\mathbf{x}) = \boldsymbol{0}$ for two-dimensional $\mathbf{f}$ and $\mathbf{x}$.

    **(b)** ⌨ Use {ref}`function-levenberg` to find an intersection point near $u=1$, $v=0.1$.

    **(d)** ⌨ Use {ref}`function-levenberg` to find an intersection point near $u=0.1$, $v=1$.

    ````{only} solutions
    ``` matlab
    f = @(x) [x(1)*log(x(1))+x(2)*log(x(2))+0.3;x(1)^4+x(2)^2-1];
    levenberg(f,[1;.1])
    levenberg(f,[.1;1])
    ```
    ````

2. ⌨ (variation of [earlier problem](problem-orbitintersect)) Two elliptical orbits $(x_1(s),y_1(s))$ and $(x_2(t),y_2(t))$ are described by the equations
  
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

    **(a)** ✍ Write out a $2\times 2$ nonlinear system of equations that describes an intersection of these orbits. (Note: An intersection is not the same as a collision---they don't have to occupy the same point at the same time.)

    **(b)** ⌨ Use {ref}`function-levenberg` to find all of the unique intersections.

    ````{only} solutions

    ``` matlab  
    x1 = @(t) -5+10*cos(t);   y1 = @(t) 6*sin(t);
    clf, fplot(x1,y1,[0 2*pi])
    x2 = @(t) 8*cos(t);   y2 = @(t) 1+12*sin(t);
    hold on, fplot(x2,y2,[0 2*pi])


    f = @(x) [-5+10*cos(x(1)) - 8*cos(x(2)); 6*sin(x(1))-1-12*sin(x(2))];
    x = levenberg(f,[-8;9]);  
    t = x(:,end);
    pos1 = [x1(t(1));y1(t(1))]
    x = levenberg(f,[8;-9]);  
    t = x(:,end);
    pos2 = [x1(t(1));y1(t(1))]
    ```
    ````

3. ⌨  (variation of [earlier problem](problem-ellipsemin)) Suppose one wants to find the points on the ellipsoid $x^2/25 + y^2/16 + z^2/9 = 1$ that are closest to and farthest from the point $(5,4,3)$. The method of Lagrange multipliers implies that any such point satisfies
  
    ```{math}
    \begin{split}
        x-5 &= \frac{\lambda x}{25} \\
        y-4 &= \frac{\lambda y}{16} \\
        z-3 &= \frac{\lambda z}{9} \\
        1 &=  \frac{1}{25}x^2 + \frac{1}{16}y^2 + \frac{1}{9}z^2
    \end{split}
    ```

    for an unknown value of $\lambda$.
  
    **(a)** Write out this system in the form $\mathbf{f}(\mathbf{u}) = \boldsymbol{0}$. (Note that the system has four variables to go with the four equations.)

    **(b)** Use {ref}`function-levenberg` with different initial guesses to find the two roots of this system. Which is the closest point to $(5,4,3)$ and which is the farthest?

    ````{only} solutions

    ``` matlab
    f = @(u) [ u(1)-5-u(1)*u(4)/25; ...
    u(2)-4-u(2)*u(4)/16; ...
    u(3)-3-u(3)*u(4)/9; ...
    u(1)^2/25 + u(2)^2/16 + u(3)^2/9 - 1 ];
    x_far = levenberg(f,[0;0;0;50])
    x_close = levenberg(f,[0;0;0;-50])
   ```
   ````

4. ✍ The Broyden update formula {eq}`broyden` is just one instance of so-called rank-one updating. Verify the  *Sherman--Morrison formula*,
  
    ```{math}
    (\mathbf{A}+\mathbf{u}\mathbf{v}^T)^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\frac{\mathbf{u}\mathbf{v}^T}{1+\mathbf{v}^T\mathbf{A}^{-1}\mathbf{u}}\mathbf{A}^{-1},
    ```

    valid whenever $\mathbf{A}$ is invertible and the denominator above is nonzero. (Hint: Show that $\mathbf{A}+\mathbf{u}\mathbf{v}^T$ times the matrix above simplifies to the identity matrix.)

5. ✍ Derive equation {eq}`nlsgradient`.
  
6. ⌨ (see also an [earlier problem](problem-newtonsysbyhand)) Suppose that

    ```{math}
    \mathbf{f}(\mathbf{x}) =
    \begin{bmatrix}
      x_1x_2+x_2^2-1 \\[1mm] x_1x_2^3 + x_1^2x_2^2 + 1
    \end{bmatrix}.
    ```

    Let $\mathbf{x}_1=[-2,1]^T$ and let $\mathbf{A}_1=\mathbf{J}(\mathbf{x}_1)$ be the exact Jacobian.
  
    **(a)** Solve {eq}`levenberg` for $\mathbf{s}_1$ with $\lambda=0$; this is the "pure" Newton step. Show numerically that $\|\mathbf{f}(\mathbf{x}_1+\mathbf{s}_1)\| > \|\mathbf{f}(\mathbf{x}_1)\|$. (Thus, the Newton step made us go to a point seemingly farther from a root than where we started.)

    **(b)** Now repeat part (a) with $\lambda=0.01j$ for $j=1,2,3,\ldots$. What is the smallest value of $j$ such that $\|\mathbf{f}(\mathbf{x}_1+\mathbf{s}_1)\| < \|\mathbf{f}(\mathbf{x}_1)\|$?

    ````{only} solutions

    ``` matlab
    f = @(x) [ x(1)*x(2)+x(2)^2-1 ; x(1)*x(2)^3 + x(1)^2*x(2)^2 + 1 ];
    A1=[1 0;-3 2];
    x1=[-2;1];
    f1 = f(x);
    lambda = 0;
    while true
        s = -(A1'*A1+lambda*eye(2)) \ (A1'*f1)
        if norm(f(x1+s)) < norm(f1), break, end
        lambda = lambda + 0.01
    end
    %% Code ends with lambda = 0.05.
    ```
    ````

7. ✍ Show that equation {eq}`levenberg` is equivalent to the linear least squares problem
  
    ```{math}
    \min_{\mathbf{v}} \Bigl(  \bigl\|\mathbf{A}_k\mathbf{v} + \mathbf{f}_k\bigr\|_2^2 +
    \lambda^2 \bigl\| \mathbf{v} \bigr\|_2^2 \Bigr).
    ```

    (Hint: Express the minimized quantity using block matrix notation.)

    Thus another interpretation of Levenberg's method is that it is the Newton step plus a penalty, weighted by $\lambda$, for taking large steps.
