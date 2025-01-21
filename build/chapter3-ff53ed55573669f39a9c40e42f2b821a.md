---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 3

Python implementations

## Functions 

(function-lsnormal-python)=
``````{dropdown} Solution of least squares by the normal equations
:open:
```{literalinclude} fncbook/fncbook/chapter03.py
:filename: lsnormal.py
:language: python
:start-line: 6
:end-line: 18
:linenos: true
```
:::{admonition} About the code
:class: dropdown
`cholesky` is imported from `scipy.linalg`.
:::
``````

(function-lsqrfact-python)=
``````{dropdown} Solution of least squares by QR factorization
:open:
```{literalinclude} fncbook/fncbook/chapter03.py
:filename: lsqrfact.py
:start-line: 20
:end-line: 30
:linenos: true
```
``````

(function-qrfact-python)=
``````{dropdown} QR factorization by Householder reflections
:open:
```{literalinclude} fncbook/fncbook/chapter03.py
:filename: qrfact.py
:start-line: 32
:end-line: 52
:linenos: true
```
``````

## Examples

```{code-cell} ipython3
:tags: remove-cell
exec(open("FNC_init.py").read())
```

### 3.1 @section-leastsq-fitting
(demo-fitting-tempinterp-python)=
``````{dropdown} @demo-fitting-tempinterp
:open:
Here are 5-year averages of the worldwide temperature anomaly as compared to the 1951–1980 average (source: NASA).

```{code-cell}
year = arange(1955,2005,5)
y = array([ -0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
    0.1180, 0.2100, 0.3320, 0.3340, 0.4560 ])

fig, ax = subplots()
ax.scatter(year, y, color="k", label="data")
xlabel("year")
ylabel("anomaly (degrees C)")
title("World temperature anomaly");
```

A polynomial interpolant can be used to fit the data. Here we build one using a Vandermonde matrix. First, though, we express time as decades since 1950, as it improves the condition number of the matrix.

```{code-cell}
t = (year - 1950) / 10
V = vander(t)
c = linalg.solve(V, y)
print(c)
```

```{index} Python; plotting functions
```

The coefficients in vector `c` are used to create a polynomial. Then we create a function that evaluates the polynomial after changing the time variable as we did for the Vandermonde matrix.

```{code-cell}
p = poly1d(c)    # convert to a polynomial
tt = linspace(1955, 2000, 500)
ax.plot(tt, p((tt - 1950) / 10), label="interpolant")
ax.legend();
fig
```

As you can see, the interpolant does represent the data, in a sense. However it's a crazy-looking curve for the application. Trying too hard to reproduce all the data exactly is known as _overfitting_.

``````

(demo-fitting-tempfit-python)=
``````{dropdown} @demo-fitting-tempfit
:open:
Here are the 5-year temperature averages again.

```{code-cell}
year = arange(1955, 2005, 5)
y = array([-0.0480, -0.0180, -0.0360, -0.0120, -0.0040,
    0.1180, 0.2100, 0.3320, 0.3340, 0.4560])
t = (year - 1950) / 10
```


The standard best-fit line results from using a linear polynomial that meets the least-squares criterion.
```{tip}
:class: dropdown
Backslash solves overdetermined linear systems in a least-squares sense.
```

```{code-cell}
V = array([ [t[i], 1] for i in range(t.size) ])    # Vandermonde-ish matrix
print(V.shape)
```

```{index} ! Python; lstsq
```

```{code-cell}
from numpy.linalg import lstsq
c, res, rank, sv = lstsq(V, y)
p = poly1d(c)
f = lambda year: p((year - 1950) / 10)
```
```

```{code-cell}
fig, ax = subplots()
ax.scatter(year, y, color="k", label="data")
yr = linspace(1955, 2000, 500)
ax.plot(yr, f(yr), label="linear fit")

xlabel("year")
ylabel("anomaly (degrees C)")
title("World temperature anomaly");
ax.legend();
```

If we use a global cubic polynomial, the points are fit more closely.

```{code-cell}
V = array([ [t[i]**3,t[i]**2,t[i],1] for i in range(t.size) ])    # Vandermonde-ish matrix
print(V.shape)
```

Now we solve the new least-squares problem to redefine the fitting polynomial.
```{tip}
:class: dropdown
The definition of `f` above is in terms of `c`. When `c` is changed, `f` is updated with it.
```

```{code-cell}
c, res, rank, sv = lstsq(V, y, rcond=None)
yr = linspace(1955, 2000, 500)
ax.plot(yr, f(yr), label="cubic fit")
fig
```

If we were to continue increasing the degree of the polynomial, the residual at the data points would get smaller, but overfitting would increase.
``````

(demo-fitting-pirate-python)=
``````{dropdown} @demo-fitting-pirate
:open:
```{code-cell}
a = array([1 / (k+1)**2 for k in range(100)])
s = cumsum(a)        # cumulative summation
p = sqrt(6*s)

plot(range(100), p, "o")
xlabel("$k$") 
ylabel("$p_k$") 
title("Sequence convergence");
```

This graph suggests that maybe $p_k\to \pi$, but it's far from clear how close the sequence gets. It's more informative to plot the sequence of errors, $\epsilon_k= |\pi-p_k|$. By plotting the error sequence on a log-log scale, we can see a nearly linear relationship.

```{code-cell}
ep = abs(pi - p)    # error sequence
loglog(range(100), ep, "o")
xlabel("$k$") 
ylabel("error") 
title("Sequence convergence");  
```

The straight line on the log-log scale suggests a power-law relationship where $\epsilon_k\approx a k^b$, or $\log \epsilon_k \approx b (\log k) + \log a$.

```{code-cell}
V = array([ [1, log(k+1)] for k in range(100) ])     # fitting matrix
c = lstsq(V, log(ep), rcond=None)[0]           # coefficients of linear fit
print(c)
```

In terms of the parameters $a$ and $b$ used above, we have

```{code-cell}
a, b = exp(c[0]), c[1]
print(f"b: {b:.3f}")
```

It's tempting to conjecture that the slope $b\to -1$ asymptotically. Here is how the numerical fit compares to the original convergence curve.

```{code-cell}
loglog(range(100), ep, "o", label="sequence")
k = arange(1,100)
plot(k, a*k**b, "--", label="power fit")
xlabel("$k$");  ylabel("error"); 
legend(); title("Sequence convergence");
```
``````

### 3.2 @section-leastsq-normaleqns

(demo-normaleqns-instab-python)=
``````{dropdown} @demo-normaleqns-instab
:open:

Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly dependent, we should find that the following matrix is somewhat ill-conditioned.

```{code-cell}
from numpy.linalg import cond
t = linspace(0, 3, 400)
A = array([ [sin(t)**2, cos((1+1e-7)*t)**2, 1] for t in t ])
kappa = cond(A)
print(f"cond(A) is {kappa:.3e}")
```

Now we set up an artificial linear least-squares problem with a known exact solution that actually makes the residual zero.

```{code-cell}
x = array([1, 2, 1])
b = A @ x
```

Using backslash to find the least-squares solution, we get a relative error that is well below $\kappa$ times machine epsilon.

```{code-cell}
from numpy.linalg import lstsq
x_BS = lstsq(A, b, rcond=None)[0]
print(f"observed error: {norm(x_BS - x) / norm(x):.3e}")
print(f"conditioning bound: {kappa * finfo(float).eps:.3e}")
```

If we formulate and solve via the normal equations, we get a much larger relative error. With $\kappa^2\approx 10^{14}$, we may not be left with more than about 2 accurate digits.

```{code-cell}
N = A.T @ A
x_NE = linalg.solve(N, A.T @ b)
relative_err = norm(x_NE - x) / norm(x)
print(f"observed error: {relative_err:.3e}")
print(f"accurate digits: {-log10(relative_err):.2f}")
```
``````

### 3.3 @section-leastsq-qr
(demo-qr-qrfact-python)=
``````{dropdown} @demo-qr-qrfact
:open:

MATLAB provides access to both the thin and full forms of the QR factorization.

```{code-cell}
A = 1.0 + floor(9 * random.rand(6,4))
A.shape
```

Here is the full form:

```{code-cell}
from numpy.linalg import qr
Q, R = qr(A, "complete")
print(f"size of Q is {Q.shape}")
print("R:")
print(R)
```

We can test that $\mathbf{Q}$ is an orthogonal matrix:

```{code-cell}
print(f"norm of (Q^T Q - I) is {norm(Q.T @ Q - eye(6)):.3e}")
```

The default for `qr`, and the one you usually want, is the thin form.

```{code-cell}
Q_hat, R_hat = qr(A)
print(f"size of Q_hat is {Q_hat.shape}")
print("R_hat:")
print(R_hat)
```

Now $\hat{\mathbf{Q}}$ cannot be an orthogonal matrix, because it is not square, but it is still ONC. Mathematically, $\hat{\mathbf{Q}}^T \hat{\mathbf{Q}}$ is a $4\times 4$ identity matrix.

```{code-cell}
print(f"norm of (Q_hat^T Q_hat - I) is {norm(Q_hat.T @ Q_hat - eye(4)):.3e}")
```
``````

(demo-qr-stable-python)=
``````{dropdown} @demo-qr-stable
:open:
We'll repeat the experiment of {numref}`Demo {number} <demo-normaleqns-instab>`, which exposed instability in the normal equations. 

```{code-cell}
t = linspace(0, 3, 400)
A = array([ [sin(t)**2, cos((1+1e-7)*t)**2, 1] for t in t ])
x = array([1, 2, 1])
b = A @ x
```

The error in the solution by {numref}`Function {number} <function-lsqrfact>` is similar to the bound predicted by the condition number.

```{code-cell}
print(f"observed error: {norm(FNC.lsqrfact(A, b) - x) / norm(x):.3e}")
print(f"conditioning bound: {cond(A) * finfo(float).eps:.3e}")
```
``````

### 3.4 @section-leastsq-house
(demo-house-qr-python)=
``````{dropdown} @demo-house-qr
:open:

We will use Householder reflections to produce a QR factorization of a matrix.

```{code-cell}
A = 1.0 + floor(9 * random.rand(6,4))
m, n = A.shape
print(A)
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 

```{code-cell}
z = A[:, 0]
v = z - norm(z) * hstack([1, zeros(m-1)])
P_1 = eye(m) - (2 / dot(v, v)) * outer(v, v)   # reflector
```

We check that this reflector introduces zeros as it should:

```{code-cell}
print(P_1 @ z)
```

Now we replace $\mathbf{A}$ by $\mathbf{P}_1\mathbf{A}$.

```{code-cell}
A = P_1 @ A
print(A)
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A[1:, 1]
v = z - norm(z) * hstack([1, zeros(m-2)])
P_2 = eye(m-1) - (2 / dot(v, v)) * outer(v, v) 
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A[1:, 1:] = P_2 @ A[1:, 1:]
print(A)
```

We need to iterate the process for the last two columns.

```{code-cell}
for j in [2, 3]:
    z = A[j:, j]
    v = z - norm(z) * hstack([1, zeros(m-j-1)])
    P = eye(m-j) - (2 / dot(v, v)) * outer(v, v)
    A[j:, j:] = P @ A[j:, j:]
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{index} Python; triu
```

```{code-cell}
R = triu(A)
print(R)
```
``````