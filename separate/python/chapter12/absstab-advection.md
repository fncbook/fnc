---
numbering:
  enumerator: 12.3.%s
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
```{code-cell}
:tags: [remove-cell]
from numpy import *
from scipy import linalg
from scipy.linalg import norm
from matplotlib.pyplot import *
from prettytable import PrettyTable
from timeit import default_timer as timer
import sys
sys.path.append('fncbook/')
import fncbook as FNC

# This (optional) block is for improving the display of plots.
# from IPython.display import set_matplotlib_formats
# set_matplotlib_formats("svg","pdf")
# %config InlineBackend.figure_format = 'svg'
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

(section-advection-absstab)=

# Absolute stability

The CFL criterion gives a necessary condition for convergence. It suggests, but cannot confirm, that a step size of $O(h)$ may be adequate in the advection equation. More details emerge when we adopt the semidiscretization point of view.

```{index} method of lines
```

Let the advection equation

```{math}
:label: advectioncc
u_t + c u_x = 0
```

be posed for $x \in [0,1]$ and subjected to periodic end conditions. If we use the central-difference matrix $\mathbf{D}_x$ defined in {eq}`trafficdiffmat` to discretize the space derivative, we get the ODE system

$$
  \mathbf{u}' = -c \mathbf{D}_x \mathbf{u}.
$$

To apply an IVP solver, we need to compare the stability region of the solver with the eigenvalues of $-c \mathbf{D}_x$, as in {numref}`section-diffusion-absstab`. You can verify (see @problem-absstab-d1eigs) that for $m$ points in $[0,1)$, these are

```{math}
:label: D1eigs
  \lambda_k = - i\, c m \sin \left( \frac{2\pi k}{m} \right), \qquad k = 0,\ldots,m-1.
```

Two things stand out about these eigenvalues: they are purely imaginary, which is consistent with conservation of magnitude, and they extend no farther than $O(m)=O(h^{-1})$ away from the origin. These characteristics suggest how to analyze the use of different time-stepping methods by referring to stability regions.

::::{prf:example} Eigenvalues for advection
:label: demo-absstab-advection

For $c=1$ we get purely imaginary eigenvalues.

```{code-cell}
:tags: hide-input
from scipy.linalg import eigvals
x, Dx, Dxx = FNC.diffper(40, [0, 1])
lamb = eigvals(Dx)

plot(real(lamb), imag(lamb), "o")
axis([-40, 40, -40, 40])
axis("equal")
title("Eigenvalues for pure advection");
```

Let's choose a time step of $\tau=0.1$ and compare to the stability regions of the Euler and backward Euler time steppers (shown as shaded regions):

```{code-cell}
:tags: hide-input
zc = exp(2j * pi * arange(361) / 360)
# points on |z|=1

z = zc - 1    # shift circle left by 1
fill(real(z), imag(z), color=(0.8, 0.8, 1))
plot(real(0.1 * lamb), imag(0.1 * lamb), "o")
axis([-5, 5, -5, 5]),  axis("equal")
title("Euler");
```

In the Euler case it's clear that *no* real value of $\tau>0$ is going to make $\zeta$ values fit within the stability region. Any method whose stability region includes none of the imaginary axis is an unsuitable choice for advection.

```{code-cell}
:tags: hide-input
z = zc + 1    # shift circle right by 1
fill([-6, 6, 6, -6], [-6, -6, 6, 6], color=(0.8, 0.8, 1))
fill(real(z), imag(z), color="w")
plot(real(0.1 * lamb), imag(0.1 * lamb), "o")
axis([-5, 5, -5, 5])
axis("equal")
title("Backward Euler");
```

The A-stable backward Euler time stepping tells the exact opposite story: it will be absolutely stable for any choice of the time step $\tau$.

::::

Many PDEs that conserve quantities will have imaginary eigenvalues, causing Euler and some other IVP methods to fail regardless of step size. Diffusion problems, in which the eigenvalues are negative and real, are compatible with a wider range of integrators, though possibly with onerous step size requirements due to stiffness.

The location of eigenvalues near $\pm ic/h$ also confirms what the CFL condition was suggesting. In order to use RK4, for example, whose stability region intersects the imaginary axis at around $\pm 2.8i$, the time step stability restriction is $\tau c/h \le 2.8$, or $\tau=O(h)$. This is much more favorable than for diffusion, whose eigenvalues were as large as $O(h^{-2})$.

```{important}
Explicit IVP methods are much more attractive for advection problems than for diffusion.
```

## Advection–diffusion equations

```{index} ! advection-diffusion equation
```

The traffic flow equation {eq}`trafficpde` combines a nonlinear advection with a diffusion term. The simplest linear problem with the same feature is the **advection–diffusion equation**

```{math}
:label: eq-advectiondiffusion
u_t+c u_x=\epsilon u_{xx}.
```

The parameter $\epsilon$ controls the relative strength between the two mechanisms, and the eigenvalues accordingly vary between the purely imaginary ones of advection and the negative real ones of diffusion.

::::{prf:example} Eigenvalues for advection–diffusion
:label: demo-absstab-advdiff

The eigenvalues of advection-diffusion are near-imaginary for $\epsilon\approx 0$ and get closer to the negative real axis as $\epsilon$ increases.

```{code-cell}
:tags: hide-input
x, Dx, Dxx = FNC.diffper(40, [0, 1])
for ep in [0.001, 0.01, 0.05]:
    lamb = eigvals(-Dx + ep * Dxx)
    plot(real(lamb), imag(lamb), "o", label=f"epsilon={ep:.1g}")
axis("equal")
legend()
title("Eigenvalues for advection-diffusion");
```

If we use a time stepping method with fixed step size $\tau$ on this system, then for any given $\epsilon$ and $m$, there is a maximum value $\hat{\tau}$ of $\tau$ such that $\hat{\tau} \lambda$ fits inside the stability region of the method for all eigenvalues $\lambda$. Then the stability restriction for the fully discrete method is $\tau \le \hat{\tau}$. (See @problem-absstab-advdiff.)

::::

In a nonlinear problem, the eigenvalues come from the linearization about an exact solution, as in {numref}`section-diffusion-stiffness`.

## Boundary effects

Boundary conditions can have a dramatic effect on the eigenvalues of the semidiscretization. For instance, @demo-upwind-direction solves linear advection $u_t=u_x$ on $[0,1]$ with the homogeneous inflow condition $u(0,t)=0$. Exclusion of the boundary node from the semidiscretization $\mathbf{u}$ to get the interior vector $\mathbf{v}$ is equivalent to

$$
\mathbf{v} = \mathbf{E} \mathbf{u}, \quad   \mathbf{u} = \begin{bmatrix}  \mathbf{v} \\ 0 \end{bmatrix} = \mathbf{E}^T \mathbf{v},
$$

where $\mathbf{E}$ is the $(m+1)\times (m+1)$ identity with the last row deleted. The ODE on the interior nodes is

$$
\frac{d\mathbf{v}}{dt} = \mathbf{E} \left( \mathbf{D}_x \mathbf{u} \right) = \mathbf{E} \mathbf{D}_x \mathbf{E}^T \mathbf{v}.
$$

As a result, we conclude that $\mathbf{A} = \mathbf{E} \mathbf{D}_x \mathbf{E}^T$ is the appropriate matrix for determining the eigenvalues of the semidiscretization. More simply, we can simply delete the last row and last column from $\mathbf{D}_x$.

::::{prf:example} Eigenvalues for an inflow boundary
:label: demo-absstab-inflow

Deleting the last row and column places all the eigenvalues of the discretization into the left half of the complex plane. 

```{code-cell}
from scipy.linalg import eigvals
x, Dx, Dxx = FNC.diffcheb(40, [0, 1])
A = -Dx[1:, 1:]  # leave out first row and column

lamb = eigvals(A)
plot(real(lamb), imag(lamb), "o")
xlim(-300, 100),  axis("equal"),  grid(True)
title("Eigenvalues of advection with zero inflow");
```

Note that the rightmost eigenvalues have real part at most

```{code-cell}
print(f"rightmost extent of eigenvalues: {max(real(lamb)):.3g}")
```

Consequently all solutions decay exponentially to zero as $t\to\infty$. This matches our observation of the solution: eventually, everything flows out of the domain.


::::

## Exercises

``````{exercise}
:label: problem-absstab-d1eigs
✍ (Similar to @problem-absstab-d2eigs.) Let $\mathbf{D}_{x}$ be $m\times m$ and given by {eq}`cflcentral` for periodic end conditions. For any integer $k \in \{0,\ldots,m-1\}$, define $\omega = \exp(2ik\pi/m)$, and let $\mathbf{v}$ be the vector whose components are $v_j = \omega^j$ for $j=0,\ldots,m-1$.

**(a)** Show that $\omega^m = 1$. 

**(b)** Let $\mathbf{v}' = \mathbf{D}_x \mathbf{v}$. Show that for $j=1,\ldots,m-2$,

$$
v_j' = \frac{1}{2h} \omega^{j} \left( \omega - \omega^{-1} \right). 
$$

**(c)** Show that the result of part (b) holds for $j=0$ and $j=m-1$ as well.

**(d)** Explain why the above results prove that $\mathbf{v}$ is an eigenvector of $\mathbf{D}_x$ with associated eigenvalue

$$
\lambda =  i\, m  \sin \left( \frac{2k\pi}{m} \right).
$$
``````

``````{exercise}
:label: problem-absstab-advdiff
⌨ In @demo-absstab-advdiff we saw that the eigenvalues of the semidiscretization of @advectioncc for periodic end conditions lie in the left half of the complex plane. Suppose we want to apply the Euler time stepping formula. For a given eigenvalue $\lambda$, there is a value of $\tau$ such that $\zeta=\tau \lambda$ lies on the boudnary of the stability region.

**(a)** ✍ Show that if $\lambda = x + iy$ and $|\zeta + 1|^2 = 1$, then  
```{math}
:label: eq-absstab-eulerperiodic
\tau = -\frac{2x}{x^2 + y^2}.
```
**(b)** ⌨ Use $c=1,$ $\epsilon=0.001,$ and $m=100$ to find the eigenvalues. For each eigenvalue, use @eq-absstab-eulerperiodic to find the value of $\tau$ that scales it to the stability region. Then find the minimum value of $\tau$ over all the eigenvalues. (This is the maximum stable time step $\hat{\tau}$ for Euler.)
``````

``````{exercise}
:label: problem-absstab-inflow
⌨ (Continuation of @problem-absstab-advdiff.) In @demo-absstab-inflow we found that the eigenvalues for @advectioncc change a great deal if an inflow boundary condition is applied. 

For $c=1$ and $m=40,$ use the reasoning in @problem-absstab-advdiff to find the maximum stable time step $\hat{\tau}$ for Euler.
``````

``````{exercise}
:label: problem-absstab-outflow
⌨ Modify @demo-absstab-inflow so that it produces the eigenvalues of the problem $u_t+u_x=0$ with an outflow condition $u(1,t)=0$. What is the behavior of solutions as $t\to\infty$?
``````
