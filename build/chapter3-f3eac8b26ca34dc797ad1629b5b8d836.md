---
kernelspec:
  display_name: MATLAB
  language: matlab
  name: jupyter_matlab_kernel
numbering:
  headings: false
---
# Chapter 3

MATLAB implementations

## Functions

(function-lsnormal-matlab)=
``````{dropdown} Solution of least squares by the normal equations
:open:
```{literalinclude} ../matlab/fnc/lsnormal.m
:language: matlab
:linenos: true
```
``````

(function-lsqrfact-matlab)=
``````{dropdown} Solution of least squares by QR factorization
:open:
```{literalinclude} ../matlab/fnc/lsqrfact.m
:language: matlab
:linenos: true
```
``````

(function-qrfact-matlab)=
``````{dropdown} QR factorization by Householder reflections
:open:
```{literalinclude} ../matlab/fnc/qrfact.m
:language: matlab
:linenos: true
```
``````

## Examples

```{code-cell}
:tags: [remove-cell]
addpath /Users/driscoll/Documents/GitHub/fnc/matlab/fnc
addpath /Users/driscoll/Documents/GitHub/fnc/matlab
FNC_init
```

### 3.1 @section-leastsq-fitting
(demo-fitting-tempinterp-matlab)=
``````{dropdown} @demo-fitting-tempinterp
Here are 5-year averages of the worldwide temperature anomaly as compared to the 1951–1980 average (source: NASA).
    
```{code-cell} matlab
t = (1955:5:2000)';
y = [ -0.0480; -0.0180; -0.0360; -0.0120; -0.0040;
    0.1180; 0.2100; 0.3320; 0.3340; 0.4560 ];
scatter(t, y), axis tight
xlabel('year')
ylabel(('anomaly ({\circ}C)'));
```

A polynomial interpolant can be used to fit the data. Here we build one using a Vandermonde matrix. First, though, we express time as decades since 1950, as it improves the condition number of the matrix.

```{code-cell} matlab
t = (t - 1950) / 10;  
n = length(t);
V = ones(n, 1);    % t^0
for j = 1:n-1
    V(:, j+1) = t .* V(:,j);
end
c = V \ y;    % solve for coefficients
```

We created the Vandermonde matrix columns in increasing-degree order. Thus, the coefficients in `c` also follow that ordering, which is the opposite of what MATLAB uses. We need to flip the coefficients before using them in `polyval`.

```{index} ! MATLAB; fplot, MATLAB; polyval
```

```{code-cell} matlab
p = @(year) polyval(c(end:-1:1), (year - 1950) / 10);
hold on
fplot(p, [1955, 2000])    % plot the interpolating function
```
``````

(demo-fitting-tempfit-matlab)=
``````{dropdown} @demo-fitting-tempfit
Here are the 5-year temperature averages again.

```{code-cell}
year = (1955:5:2000)';
y = [ -0.0480; -0.0180; -0.0360; -0.0120; -0.0040;
    0.1180; 0.2100; 0.3320; 0.3340; 0.4560 ];
```

```{index} ! MATLAB; \\
```

::::{grid} 1 1 2 2
The standard best-fit line results from using a linear polynomial that meets the least-squares criterion.
:::{card}
Backslash solves overdetermined linear systems in a least-squares sense.
:::
::::

```{code-cell}
t = (year - 1955) / 10;    % better matrix conditioning later
V = [ t.^0 t ];            % Vandermonde-ish matrix
size(V)
```

```{code-cell}
c = V \ y;
f = @(year) polyval(c(end:-1:1), (year - 1955) / 10);
```

```{code-cell}
clf
scatter(year, y), axis tight
xlabel('year'), ylabel('anomaly ({\circ}C)')
hold on
fplot(f, [1955, 2000]) 
```

If we use a global cubic polynomial, the points are fit more closely.

```{code-cell}
V = [t.^0, t.^1, t.^2, t.^3];    % Vandermonde-ish matrix  
size(V)
```

::::{grid} 1 1 2 2
Now we solve the new least-squares problem to redefine the fitting polynomial.
:::{card}
The definition of `f` above is in terms of `c`. When `c` is changed, then `f` has to be redefined.
:::
::::

```{code-cell}
c = V \ y;
f = @(year) polyval(c(end:-1:1), (year - 1955) / 10);
fplot(f, [1955, 2000]) 
legend('data', 'linear', 'cubic', 'Location', 'northwest');
```

If we were to continue increasing the degree of the polynomial, the residual at the data points would get smaller, but overfitting would increase.
``````

(demo-fitting-pirate-matlab)=
``````{dropdown} @demo-fitting-pirate
```{code-cell}
k = (1:100)';
a = 1./k.^2;      % sequence
s = cumsum(a);    % cumulative summation
p = sqrt(6*s);
clf
plot(k, p, 'o-')
xlabel('k'), ylabel('p_k')
title(('Sequence converging to \pi'));
```

This graph suggests that maybe $p_k\to \pi$, but it's far from clear how close the sequence gets. It's more informative to plot the sequence of errors, $\epsilon_k= |\pi-p_k|$. By plotting the error sequence on a log-log scale, we can see a nearly linear relationship.

```{code-cell}
ep = abs(pi - p);    % error sequence
loglog(k, ep, 'o')
title('Convergence')
xlabel('k'), ylabel('|p_k - \pi|'), axis tight    
```

The straight line on the log-log scale suggests a power-law relationship where $\epsilon_k\approx a k^b$, or $\log \epsilon_k \approx b (\log k) + \log a$.

```{code-cell}
V = [ k.^0, log(k) ];    % fitting matrix
c = V \ log(ep)          % coefficients of linear fit
```

In terms of the parameters $a$ and $b$ used above, we have

```{code-cell}
a = exp(c(1)),  b = c(2)
```

It's tempting to conjecture that the slope $b\to -1$ asymptotically. Here is how the numerical fit compares to the original convergence curve.

```{code-cell}
hold on
loglog(k, a * k.^b)
legend('sequence', 'power-law fit');
```
``````

### 3.2 @section-leastsq-normaleqns
(demo-normaleqns-instab-matlab)=
``````{dropdown} @demo-normaleqns-instab

::::{grid} 1 1 2 2
Because the functions $\sin^2(t)$, $\cos^2(t)$, and $1$ are linearly dependent, we should find that the following matrix is somewhat ill-conditioned.
:::{card}
The local variable scoping rule for loops applies to comprehensions as well.
:::
::::

```{code-cell}
t = linspace(0, 3, 400)';
A = [ sin(t).^2, cos((1+1e-7)*t).^2, t.^0 ];
kappa = cond(A)
```

Now we set up an artificial linear least-squares problem with a known exact solution that actually makes the residual zero.

```{code-cell}
x = [1; 2; 1];
b = A * x;
```

Using backslash to find the least-squares solution, we get a relative error that is well below $\kappa$ times machine epsilon.

```{code-cell}
x_BS = A \ b;
observed_err = norm(x_BS - x) / norm(x)
max_err = kappa * eps
```

If we formulate and solve via the normal equations, we get a much larger relative error. With $\kappa^2\approx 10^{14}$, we may not be left with more than about 2 accurate digits.

```{code-cell}
N = A'*A;
x_NE = N\(A'*b);
observed_err = norm(x_NE - x) / norm(x)
digits = -log10(observed_err)
```
``````

### 3.3 @section-leastsq-qr
(demo-qr-qrfact-matlab)=
``````{dropdown} @demo-qr-qrfact

MATLAB provides access to both the thin and full forms of the QR factorization.

```{code-cell}
A = magic(5);
A = A(:, 1:4);
[m, n] = size(A)
```

Here is the full form:

```{code-cell}
[Q, R] = qr(A);
szQ = size(Q), szR = size(R)
```

We can test that $\mathbf{Q}$ is an orthogonal matrix:

```{code-cell}
QTQ = Q' * Q
norm(QTQ - eye(m))
```

With a second input argument given to `qr`, the thin form is returned. (This is usually the one we want in practice.)

```{code-cell}
[Q_hat, R_hat] = qr(A, 0);
szQ_hat = size(Q_hat), szR_hat = size(R_hat)
```

Now $\hat{\mathbf{Q}}$ cannot be an orthogonal matrix, because it is not square, but it is still ONC. Mathematically, $\hat{\mathbf{Q}}^T \hat{\mathbf{Q}}$ is a $4\times 4$ identity matrix.

```{code-cell}
Q_hat' * Q_hat - eye(n)
```
``````

(demo-qr-stable-matlab)=
``````{dropdown} @demo-qr-stable
We'll repeat the experiment of {numref}`Demo {number} <demo-normaleqns-instab>`, which exposed instability in the normal equations. 

```{code-cell}
t = linspace(0, 3, 400)';
A = [ sin(t).^2, cos((1+1e-7)*t).^2, t.^0 ];
x = [1; 2; 1];
b = A * x;
```

The error in the solution by {numref}`Function {number} <function-lsqrfact>` is similar to the bound predicted by the condition number.

```{code-cell}
observed_error = norm(lsqrfact(A, b) - x) / norm(x)
error_bound = cond(A) * eps
```
``````

### 3.4 @section-leastsq-house

(demo-house-qr-matlab)=
``````{dropdown} @demo-house-qr

We will use Householder reflections to produce a QR factorization of a matrix.

```{code-cell}
A = magic(6);
A = A(:, 1:4);
[m, n] = size(A)
```

```{index} ! MATLAB; eye
```

Our first step is to introduce zeros below the diagonal in column 1 by using {eq}`hhvector` and {eq}`hhreflect`. 

```{code-cell}
z = A(:, 1);
v = z - norm(z) * eye(m,1);
P_1 = eye(m) - 2 / (v' * v) * (v * v'); 
```

We check that this reflector introduces zeros as it should:

```{code-cell}
P_1 * z
```

Now we replace $\mathbf{A}$ by $\mathbf{P}_1\mathbf{A}$.

```{code-cell}
A = P_1 * A
```

We are set to put zeros into column 2. We must not use row 1 in any way, lest it destroy the zeros we just introduced. So we leave it out of the next reflector.

```{code-cell}
z = A(2:m, 2);
v = z - norm(z) * eye(m-1, 1);
P_2 = eye(m-1) - 2 / (v' * v) * (v * v');
```

We now apply this reflector to rows 2 and below only.

```{code-cell}
A(2:m, 2:n) = P_2 * A(2:m, 2:n)
```

We need to iterate the process for the last two columns.

```{code-cell}
for j = 3:n
    z = A(j:m,j);
    k = m-j+1;
    v = z - norm(z) * eye(k, 1);
    P = eye(k) - 2 / (v' * v) * (v * v');
    A(j:m, j:n) = P * A(j:m, j:n);
end
```

We have now reduced the original to an upper triangular matrix using four orthogonal Householder reflections:

```{code-cell}
R = A
```
``````
