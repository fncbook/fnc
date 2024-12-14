---
kernelspec:
  display_name: Python 3
  language: python
  name: python3
numbering:
  headings: false
---
# Chapter 5

Python implementations

## Functions 

(function-hatfun-python)=
``````{dropdown} Hat function
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: hatfun.py
:start-line: 3
:end-line: 25
:language: python
:linenos: true
```
``````

(function-plinterp-python)=
``````{dropdown} Piecewise linear interpolation
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: plinterp.py
:start-line: 27
:end-line: 35
:language: python
:linenos: true
```
``````

(function-spinterp-python)=
``````{dropdown} Cubic spline interpolation
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: spinterp.py
:start-line: 37
:end-line: 95
:language: python
:linenos: true
```
``````

(function-fdweights-python)=
``````{dropdown} Fornberg's algorithm for finite difference weights
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: fdweights.py
:start-line: 97
:end-line: 138
:language: python
:linenos: true
```
``````

(function-trapezoid-python)=
``````{dropdown} Trapezoid formula for numerical integration
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: trapezoid.py
:start-line: 140
:end-line: 150
:language: python
:linenos: true
```
``````

(function-intadapt-python)=
``````{dropdown} Adaptive integration
```{literalinclude} ../python/pkg/FNC/FNC05.py
:filename: intadapt.py
:start-line: 152
:end-line: 191
:language: python
:linenos: true
```
:::{admonition} About the code
:class: dropdown
The intended way for a user to call {numref}`Function {number} <function-intadapt>` is with only `f`, `a`, `b`, and `tol` provided. We then use default values on the other parameters to compute the function values at the endpoints, the interval's midpoint, and the function value at the midpoint. Recursive calls from within the function itself will provide all of that information, since it was already calculated along the way.
:::
``````

## Examples

```{code-cell} ipython3
import FNC
from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import solve, norm
import scipy.sparse as sparse
from scipy.sparse.linalg import splu
from timeit import default_timer as timer
```

```{code-cell} ipython3
:tags: [remove-cell]
# This (optional) block is for improving the display of plots.
rcParams["figure.figsize"] = [7, 4]
rcParams["lines.linewidth"] = 2
rcParams["lines.markersize"] = 4
rcParams['animation.html'] = "jshtml"  # or try "html5"
```

```{code-cell} ipython3
f = lambda x: exp(sin(7*x))
a = 0;  b = 2;
T,t,y = FNC.trapezoid(f,a,b,40)
print("Trapezoid:",T)
```