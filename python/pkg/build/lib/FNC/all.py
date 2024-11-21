def horner(c,x):
    """
    horner(c,x)

    Evaluate a polynomial whose coefficients are given in descending order in `c`, at the point `x`, using Horner's rule.
    """

    n = len(c)
    y = c[0]
    for k in range(1,n):
        y = x*y + c[k]   

    return yfrom scipy import *

def forwardsub(L,b):
	"""
 	forwardsub(L,b)

	Solve the lower-triangular linear system with matrix L and right-hand side
	vector b.
	"""
	n = len(b)
	x = zeros(n)
	for i in range(n):
		s = L[i,:i] @ x[:i]
		x[i] = ( b[i] - s ) / L[i,i]
	return x


def backsub(U,b):
	"""
	backsub(U,b)

	Solve the upper-triangular linear system with matrix U and right-hand side
	vector b.
	"""
	n = len(b)
	x = zeros(n)
	for i in range(n-1,-1,-1):
		s = U[i,i+1:] @ x[i+1:]
		x[i] = ( b[i] - s ) / U[i,i]
	return x

def lufact(A):
	"""
	lufact(A)

	Compute the LU factorization of square matrix A, returning the factors.
	"""
	n = A.shape[0]
	L = eye(n)      # puts ones on diagonal
	U = copy(A)

	# Gaussian elimination
	for j in range(n-1):
		for i in range(j+1,n):
			L[i,j] = U[i,j] / U[j,j]   # row multiplier
			U[i,j:] = U[i,j:] - L[i,j]*U[j,j:]
	return L,triu(U)
from scipy import *
from scipy.linalg import qr

def lsnormal(A,b):
	"""
    lsnormal(A,b)

	Solve a linear least squares problem by the normal equations. Returns the
	minimizer of ||b-Ax||.
	"""

	N = A.T@A
	z = A.T@b
	R = cholesky(N)
	w = forwardsub(R.T,z)                   # solve R'z=c
	x = backsub(R,w)                        # solve Rx=z

	return x

def lsqrfact(A,b):
	"""
	lsqrfact(A,b)

	Solve a linear least squares problem by QR factorization. Returns the
	minimizer of ||b-Ax||.
	"""

	Q,R = qr(A)
	c = Q.T@b
	x = backsub(R,c)

	return x
from scipy import *
from scipy.linalg import norm,lstsq
from numpy import finfo
eps = finfo(double).eps 
import warnings

def newton(f,dfdx,x1):
	"""
	newton(f,dfdx,x1)

	Use Newton's method to find a root of `f` starting from `x1`, where `dfdx` is the
	derivative of `f`. Returns a vector of root estimates.
	"""
	# Operating parameters.
	funtol = 100*eps
	xtol = 100*eps  
	maxiter = 40

	x = zeros(maxiter)
	x[0] = x1
	y = f(x1)
	dx = Inf   # for initial pass below
	k = 0

	while (abs(dx) > xtol) and (abs(y) > funtol) and (k < maxiter):
		dydx = dfdx(x[k])
		dx = -y/dydx            # Newton step
		x[k+1] = x[k]+dx        # new estimate

		k = k+1
		y = f(x[k])

	if k==maxiter:
		warnings.warn("Maximum number of iterations reached.")

	return x[:k+1]

def secant(f,x1,x2):
	"""
	secant(f,x1,x2)

	Use the secant method to find a root of `f` starting from `x1` and `x2`. Returns a
	vector of root estimates.
	"""
    # Operating parameters.
	funtol = 100*eps
	xtol = 100*eps  
	maxiter = 40

	x = zeros(maxiter)
	x[:2] = [x1,x2]
	y1 = f(x1)
	y2 = 100
	dx = Inf   # for initial pass below
	k = 1

	while (abs(dx) > xtol) and (abs(y2) > funtol) and (k < maxiter):
		y2 = f(x[k])
		dx = -y2 * (x[k]-x[k-1]) / (y2-y1)   # secant step
		x[k+1] = x[k]+dx        # new estimate

		k = k+1
		y1 = y2    # current f-value becomes the old one next time

	if k==maxiter:		
		warnings.warn("Maximum number of iterations reached.")

	return x[:k+1]

def newtonsys(f,x1):
	"""
		newtonsys(f,x1)

	Use Newton's method to find a root of a system of equations, starting from `x1`. The
	function `f` should return both the residual vector and the Jacobian matrix. Returns
	root estimates as a matrix, one estimate per column.
	"""
	# Operating parameters.
	funtol = 1000*eps
	xtol = 1000*eps
	maxiter = 40

	x = zeros((len(x1),maxiter))
	x[:,0] = x1
	y,J = f(x1)
	dx = 10.   # for initial pass below
	k = 0

	while (norm(dx) > xtol) and (norm(y) > funtol) and (k < maxiter):
		dx = -lstsq(J,y)[0]            # Newton step
		x[:,k+1] = x[:,k] + dx

		k = k+1
		y,J = f(x[:,k])

	if k==maxiter:
		warnings.warn("Maximum number of iterations reached.")

	return x[:,:k+1]

def fdjac(f,x0,y0):
	"""
	fdjac(f,x0,y0)

	Compute a finite-difference approximation of the Jacobian matrix for `f` at `x0`,
	where `y0`=`f(x0)` is given.
	"""

	delta = sqrt(eps)   # FD step size
	m,n = (len(y0),len(x0))
	J = zeros((m,n))
	I = eye(n)
	for j in range(n):
		J[:,j] = ( f(x0 + delta*I[:,j]) - y0) / delta

	return J

def levenberg(f,x1,tol=1e-12):
	"""
	levenberg(f,x1,tol)

	Use Levenberg's quasi-Newton iteration to find a root of the system `f`, starting from
	`x1`, with `tol` as the stopping tolerance in both step size and residual norm. Returns
	root estimates as a matrix, one estimate per column.
	"""

	# Operating parameters.
	ftol = tol
	xtol = tol
	maxiter = 40

	n = len(x1)
	x = zeros((n,maxiter))
	x[:,0] = x1
	fk = f(x1)
	k = 0
	s = 10.
	Ak = fdjac(f,x[:,0],fk)   # start with FD Jacobian
	jac_is_new = True

	lam = 10
	while (norm(s) > xtol) and (norm(fk) > ftol) and (k < maxiter):
		# Compute the proposed step.
		B = Ak.T@Ak + lam*eye(n)
		z = Ak.T@fk
		s = -lstsq(B,z)[0]

		xnew = x[:,k] + s;   
		fnew = f(xnew)

		# Do we accept the result?
		if norm(fnew) < norm(fk):    # accept
			y = fnew - fk
			x[:,k+1] = xnew
			fk = fnew
			k = k+1

			lam = lam/10   # get closer to Newton
			# Broyden update of the Jacobian.
			Ak = Ak + outer( y-Ak@s, s/dot(s,s) )
			jac_is_new = False
		else:                       # don't accept
			# Get closer to steepest descent.
			lam = lam*4
			# Re-initialize the Jacobian if it's out of date.
			if not jac_is_new:
				Ak = fdjac(f,x[:,k],fk)
				jac_is_new = True

	if (norm(fk) > 1e-3):
		warnings.warn("Iteration did not find a root.")

	return x[:,:k+1]
from scipy import *
from numpy import *
from scipy.linalg import *
from numpy.linalg import *

def hatfun(x,t,k):
	"""
    hatfun(x,t,k)

	Evaluate a piecewise linear "hat" function at `x`, where `t` is a vector of
	n+1 interpolation nodes and `k` is an integer in 0:n giving the index of the node
	where the hat function equals one.
	"""

	n = len(t)-1
    # Return correct node given mathematical index k, including fictitious choices.   
	def node(k):
		if k < 0:
			return 2*t[0]-t[1]
		elif k > n: 
			return 2*t[n]-t[n-1] 
		else:
			return t[k]

	H1 = (x-node(k-1)) / (node(k)-node(k-1))   # upward slope
	H2 = (node(k+1)-x) / (node(k+1)-node(k))   # downward slope
	H = minimum(H1,H2)
	return maximum(0,H)

def plinterp(t,y):
	"""
	plinterp(t,y)

	Create a piecewise linear interpolating function for data values in `y` given at nodes
	in `t`.
	"""

	n = len(t)-1
	return lambda x: sum( y[k]*hatfun(x,t,k) for k in range(n+1) )



def spinterp(t,y):
	"""
	spinterp(t,y)

	Create a cubic not-a-knot spline interpolating function for data values in `y` given at nodes in `t`.
	"""
	n = len(t) - 1
	h = [t[i+1]-t[i] for i in range(n)]    

	# Preliminary definitions.
	Z = zeros([n,n])
	I = eye(n);  E = I[:n-1,:]
	J = eye(n) + diag(-ones(n-1),1)
	H = diag(h)

	# Left endpoint interpolation:
	AL = hstack([I,Z,Z,Z])
	vL = y[:-1]

	# Right endpoint interpolation:
	AR = hstack([I,H,H**2,H**3])
	vR = y[1:]

	# Continuity of first derivative:
	A1 = E @ hstack([Z,J,2*H,3*H**2])
	v1 = zeros(n-1)

	# Continuity of second derivative:
	A2 = E @ hstack([Z,Z,J,3*H])
	v2 = zeros(n-1)

	# Not-a-knot conditions:
	nakL = hstack([ zeros(3*n), hstack([1,-1,zeros(n-2)]) ])
	nakR = hstack([ zeros(3*n), hstack([zeros(n-2),1,-1]) ])

    # Assemble and solve the full system.
	A = vstack( [ AL, AR, A1, A2, nakL, nakR ] )
	v = hstack( [ vL, vR, v1, v2, 0, 0 ] )
	z = solve(A,v)

	# Break the coefficients into separate vectors.
	rows = arange(n)
	a = z[rows]
	b = z[n+rows];  c = z[2*n+rows];  d = z[3*n+rows]
	S = [ poly1d([d[k],c[k],b[k],a[k]]) for k in range(n) ]
	
	# This function evaluates the spline when called with a value for x.
	def evaluate(x):
		f = zeros(x.shape)
		for k in range(n): 
			# Evaluate this piece's cubic at the points inside it.
			index = (x>=t[k]) & (x<=t[k+1])
			f[index] = S[k](x[index]-t[k])
		return f

	return evaluate

def fdweights(t,m):
	"""
	fdweights(t,m)

	Return weights for the `m`th derivative of a function at zero using values at the
	nodes in vector `t`.
	"""
    # This is a compact implementation, not an efficient one.

	def weight(t,m,r,k):
		# Recursion for one weight.
		# Input:
		#   t   nodes (vector)
		#   m   order of derivative sought
		#   r   number of nodes to use from t (<= length(t))
		#   k   index of node whose weight is found

		if (m<0) or (m>r):        # undefined coeffs must be zero
			c = 0
		elif (m==0) and (r==0):  # base case of one-point interpolation
			c = 1
		else:                     # generic recursion
			if k<r:
				c = (t[r]*weight(t,m,r-1,k) -
					m*weight(t,m-1,r-1,k))/(t[r]-t[k])
			else:
				if r <= 1:
					numer = 1.0
				else:
					numer = prod(t[r-1]-t[:r-1])
				if r <= 0: 
					denom = 1.0
				else:
					denom = prod(t[r]-t[:r]) 
				beta =  numer/denom
				c = beta*(m*weight(t,m-1,r-1,r-1) - t[r-1]*weight(t,m,r-1,r-1))
		return c

	r = len(t)-1
	w = zeros(t.shape)
	return [ weight(t,m,r,k) for k in range(r+1) ] 
	#return weight(t,m,r,0)

def trapezoid(f,a,b,n):
	"""
	trapezoid(f,a,b,n)

	Apply the trapezoid integration formula for integrand `f` over interval [`a`,`b`], broken up into `n` equal pieces. Returns estimate, vector of nodes, and vector of integrand values at the nodes.
	"""
	h = (b-a)/n
	t = linspace(a,b,n+1)
	y = f(t)
	T = h * ( sum(y[1:-1]) + 0.5*(y[0] + y[-1]) )

	return T,t,y

def intadapt(f,a,b,tol):
	"""
	intadapt(f,a,b,tol)

	Do adaptive integration to estimate the integral of `f` over [`a`,`b`] to desired
	error tolerance `tol`. Returns estimate and a vector of evaluation nodes used.
	"""
    # Use error estimation and recursive bisection.
	def do_integral(a,fa,b,fb,m,fm,tol):
		# These are the two new nodes and their f-values.
		xl = (a+m)/2;  fl = f(xl)
		xr = (m+b)/2;  fr = f(xr)
		t = array([a,xl,m,xr,b])             # all 5 nodes at this level

		# Compute the trapezoid values iteratively.
		h = (b-a)
		T = array([0.,0.,0.])
		T[0] = h*(fa+fb)/2
		T[1] = T[0]/2 + (h/2)*fm
		T[2] = T[1]/2 + (h/4)*(fl+fr)

		S = (4*T[1:]-T[:-1]) / 3      # Simpson values
		E = (S[1]-S[0]) / 15           # error estimate

		if abs(E) < tol*(1+abs(S[1])):  # acceptable error?
			Q = S[1]                    # yes--done
		else:
			# Error is too large--bisect and recurse.
			QL,tL = do_integral(a,fa,m,fm,xl,fl,tol)
			QR,tR = do_integral(m,fm,b,fb,xr,fr,tol)
			Q = QL + QR
			t = hstack([tL,tR[1:]])    # merge the nodes w/o duplicate
		return Q,t

	m = (b+a)/2
	Q,t = do_integral(a,f(a),b,f(b),m,f(m),tol)
	return Q,t
from scipy import *
from numpy import *
from scipy.linalg import *
from numpy.linalg import *
from FNC04 import levenberg
import warnings

def eulerivp(dudt,tspan,u0,n):
	"""
	eulerivp(dudt,tspan,u0,n)

	Apply Euler's method to solve the IVP u'=`dudt`(u,t) over the interval `tspan` with
	u(`tspan[1]`)=`u0`, using `n` subintervals/steps. Return vectors of times and solution
	values.
	"""
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)
	u = zeros(n+1)
	u[0] = u0
	for i in range(n):
		u[i+1] = u[i] + h*dudt(t[i],u[i])
	return t,u

def eulersys(dudt,tspan,u0,n):
	"""
	eulersys(dudt,tspan,u0,n)

	Apply Euler's method to solve the vector-valued IVP u'=`dudt`(u,p,t) over the interval
	`tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
	"""
	# Time discretization.
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)

	# Initial condition and output setup.
	u = zeros([u0.size,n+1])
	u[:,0] = u0

	# The time stepping iteration.
	for i in range(n):
		u[:,i+1] = u[:,i] + h*dudt(t[i],u[:,i])

	return t,u


def ie2(dudt,tspan,u0,n):
	"""
	ie2(dudt,tspan,u0,n)

	Apply the Improved Euler method to solve the vector-valued IVP u'=`dudt`(u,p,t) over the
	interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps. Returns a vector
	of times and a vector of solution values/vectors.
	"""
	# Time discretization.
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)

	# Initialize output.
	u = zeros([len(u0),n+1])
	u[:,0] = u0

	# Time stepping.
	for i in range(n):
		uhalf = u[:,i] + h/2*dudt(t[i],u[:,i])
		u[:,i+1] = u[:,i] + h*dudt(t[i]+h/2,uhalf)

	return t,u

def rk4(dudt,tspan,u0,n):
	"""
	rk4(dudt,tspan,u0,n)

	Apply "the" Runge-Kutta 4th order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
	over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
	Return a vector of times and a vector of solution values/vectors.
	"""
	# Time discretization.
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)

	# Initialize output.
	u = zeros([len(u0),n+1])
	u[:,0] = u0

	# Time stepping.
	for i in range(n):
		k1 = h*dudt( t[i]     , u[:,i],      )
		k2 = h*dudt( t[i]+h/2 , u[:,i]+k1/2  )
		k3 = h*dudt( t[i]+h/2 , u[:,i]+k2/2  )
		k4 = h*dudt( t[i]+h   , u[:,i]+k3    )
		u[:,i+1] = u[:,i] + (k1 + 2*(k2 + k3) + k4)/6
	return t,u

def rk23(dudt,tspan,u0,tol):
	"""
	rk23(dudt,tspan,u0,tol)

	Apply adaptive embedded RK formula to solve the vector-valued IVP u'=`dudt`(u,p,t)
	over the interval `tspan` with u(`tspan[1]`)=`u0`, with error tolerance `tol`.
	Return a vector of times and a vector of solution values/vectors.
	"""
	# Initialize for the first time step.
	t = [tspan[0]]
	u = [u0];   i = 0;
	h = 0.5*tol**(1/3)
	s1 = dudt(t,u0)

	# Time stepping.
	while t[i] < tspan[-1]:
		# Detect underflow of the step size.
		if t[i]+h == t[i]:
			warnings.warn(f"Stepsize too small near t={t[i]}")
			break  # quit time stepping loop

		# New RK stages.
		s2 = dudt( t[i]+h/2,   u[i]+(h/2)*s1   )
		s3 = dudt( t[i]+3*h/4, u[i]+(3*h/4)*s2 )
		unew2 = u[i] + h*(2*s1 + 3*s2 + 4*s3)/9   # 2rd order solution
		s4 = dudt( t[i]+h , unew2 )
		err = h*(-5*s1/72 + s2/12 + s3/9 - s4/8)    # 2nd/3rd order difference
		E = norm(err,Inf)                           # error estimate
		maxerr = tol*(1 + norm(u[i],Inf))         # relative/absolute blend

		# Accept the proposed step?
		if E < maxerr:     # yes
			t.append(t[i] + h)
			u.append(unew2)
			i = i+1
			s1 = s4      # use FSAL property

		# Adjust step size.
		q = 0.8*(maxerr/E)**(1/3)       # conservative optimal step factor
		q = min(q,4)                    # limit stepsize growth
		h = min(q*h,tspan[-1]-t[i])      # don't step past the end

	# Convert outputs to arrays
	return array(t),array(u).T

def ab4(dudt,tspan,u0,n):
	"""
	ab4(dudt,tspan,u0,n)

	Apply the Adams-Bashforth 4th order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
	over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
	"""
	# Time discretization.
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)

	# Constants in the AB4 method.
	k = 4;    sigma = array([55, -59, 37, -9])/24;

	# Find starting values by RK4.
	ts,us = rk4(dudt,[a,a+(k-1)*h],u0,k-1)
	u = zeros([u0.size,n+1])
	u[:,:k] = us[:,:k]

	# Compute history of u' values, from newest to oldest.
	f = array([ dudt(t[k-j-2],u[:,k-j-2]) for j in range(k)  ])
	
	# Time stepping.
	for i in range(k-1,n):
		f = vstack([dudt(t[i],u[:,i]),f[:-1]])  # new value of du/dt
		u[:,i+1] = u[:,i] + h*dot(sigma,f)       # advance one step
	return t,u

def am2(dudt,tspan,u0,n):
	"""
	am2(dudt,tspan,u0,n)

	Apply the Adams-Moulton 2nd order method to solve the vector-valued IVP u'=`dudt`(u,p,t)
	over the interval `tspan` with u(`tspan[1]`)=`u0`, using `n` subintervals/steps.
	"""
	# Time discretization.
	a,b = tspan
	h = (b-a)/n
	t = linspace(a,b,n+1)

	# Initialize output.
	u = zeros([u0.size,n+1])
	u[:,0] = u0

	# Time stepping.
	for i in range(n):
		# Data that does not depend on the new value.
		known = u[:,i] + h/2*dudt(t[i],u[:,i])
		# Find a root for the new value.
		F = lambda z: z - h/2*dudt(t[i+1],z) - known
		unew = levenberg(F,known)
		u[:,i+1] = unew[:,-1]

	return t,u
from scipy import *
from numpy import *
from scipy.linalg import *
from numpy.linalg import *
from scipy.sparse import csc_matrix,diags

def poweriter(A,numiter):
	"""
	poweriter(A,numiter)

	Perform `numiter` power iterations with the matrix `A`, starting from a random vector,	and return a vector of eigenvalue estimates and the final eigenvector approximation.
	"""
	n = A.shape[0]
	x = randn(n)
	x = x/norm(x,Inf)
	gamma = zeros(numiter)
	for k in range(numiter):
		y = A@x
		m = argmax(abs(y))
		gamma[k] = y[m]/x[m]
		x = y/y[m]

	return gamma,x

def inviter(A,s,numiter):
	"""
	inviter(A,s,numiter)

	Perform `numiter` inverse iterations with the matrix `A` and shift `s`, starting
	from a random vector, and return a vector of eigenvalue estimates and the final
	eigenvector approximation.
	"""
	n = A.shape[0]
	x = randn(n)
	x = x/norm(x,Inf)
	gamma = zeros(numiter)
	PL,U = lu(A - s*eye(n),permute_l=True)
	for k in range(numiter):
		y = solve(U,solve(PL,x))
		m = argmax(abs(y))
		gamma[k] = x[m]/y[m] + s
		x = y/y[m]

	return gamma,x

def arnoldi(A,u,m):
	"""
	arnoldi(A,u,m)

	Perform the Arnoldi iteration for `A` starting with vector `u`, out to the Krylov
	subspace of degree `m`. Return the orthonormal basis (`m`+1 columns) and the upper
	Hessenberg `H` of size `m`+1 by `m`.
	"""
	n = u.size
	Q = zeros([n,m+1])
	H = zeros([m+1,m])
	Q[:,0] = u/norm(u)
	for j in range(m):
		# Find the new direction that extends the Krylov subspace.
		v = A@Q[:,j]
		# Remove the projections onto the previous vectors.
		for i in range(j+1):
			H[i,j] = Q[:,i]@v
			v -= H[i,j]*Q[:,i]
		# Normalize and store the new basis vector.
		H[j+1,j] = norm(v)
		Q[:,j+1] = v/H[j+1,j]

	return Q,H

def arngmres(A,b,m):
	"""
	arngmres(A,b,m)

	Do `m` iterations of GMRES for the linear system `A`*x=`b`. Return the final solution
	estimate x and a vector with the history of residual norms. (This function is for
	demo only, not practical use.)
	"""
	n = len(b)
	Q = zeros([n,m+1])
	Q[:,0] = b/norm(b)
	H = zeros([m+1,m])

	# Initial "solution" is zero.
	residual = hstack([norm(b),zeros(m)])

	for j in range(m):
		# Next step of Arnoldi iteration.
		v = A@Q[:,j]
		for i in range(j+1):
			H[i,j] = Q[:,i]@v
			v -= H[i,j]*Q[:,i]
		H[j+1,j] = norm(v)
		Q[:,j+1] = v/H[j+1,j]

		# Solve the minimum residual problem.
		r = hstack([norm(b), zeros(j+1)])
		z = lstsq(H[:j+2,:j+1],r)[0]
		x = Q[:,:j+1]@z
		residual[j+1] = norm( A@x - b )

	return x,residual

def sprandsym(n,density,**kwargs):
	if "rcond" in kwargs:
		ev = array([ kwargs["rcond"]**(i/(n-1)) for i in range(n) ])
	elif "eigvals" in kwargs:
		ev = kwargs["eigvals"]
	else:
		ev = sqrt(n)*randn(n)

	def randjr(A):
		# Random Jacobi rotation similarity transformation.
		theta = 2*pi*rand()
		c = cos(theta); s = sin(theta)
		idx = [int(k) for k in floor(n*rand(2)) ]
		R = array([[c,s],[-s,c]])
		A[idx,:] = R @ A[idx,:]
		A[:,idx] = A[:,idx] @ R.T
		return A

	targetnz = ceil(min(0.98,density)*n*n)
	A = diags(ev,0,format="lil")
	while A.nnz < targetnz:
		A = randjr(A)
	
	return csc_matrix(A)
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from scipy.linalg import *
from numpy.linalg import *

def polyinterp(t,y):
    """
    polyinterp(t,y)

    Return a callable polynomial interpolant through the points in vectors `t`,`y`. Uses
    the barycentric interpolation formula.
    """
    n = len(t)-1
    C = (t[-1]-t[0]) / 4       # scaling factor to ensure stability
    tc = t/C

    # Adding one node at a time, compute inverses of the weights.
    omega = ones(n+1)
    for m in range(n):
        d = tc[:m+1] - tc[m+1]          # vector of node differences
        omega[:m+1] = omega[:m+1]*d     # update previous
        omega[m+1] = prod( -d )         # compute the new one
    w = 1 / omega                   # go from inverses to weights

    def p(x):
        # Compute interpolant.
        z = where(x==t)[0]
        if len(z)>0:       # avoid dividing by zero
            # Apply L'Hopital's Rule exactly.
             f = y[z[0]]
        else:
            terms = w / (x - t)
            f = sum(y*terms) / sum(terms)
        return f
        
    return vectorize(p)

def triginterp(t,y):
    """
        triginterp(t,y)

    Return trigonometric interpolant for points defined by vectors `t` and `y`.
    """
    N = len(t)

    def trigcardinal(x):
        if x==0:
            tau = 1.
        elif mod(N,2)==1:      # odd
            tau = sin(N*pi*x/2) / (N*sin(pi*x/2))
        else:                # even
            tau = sin(N*pi*x/2) / (N*tan(pi*x/2))
        return tau

    def p(x):
        return sum( [y[k]*trigcardinal(x-t[k]) for k in range(N)] )
        
    return vectorize(p)

def ccint(f,n):
    """
    ccint(f,n)

    Perform Clenshaw-Curtis integration for the function `f` on `n`+1 nodes in [-1,1]. Return
    integral and a vector of the nodes used. Note: `n` must be even.
    """
    # Find Chebyshev extreme nodes.
    theta = linspace(0,pi,n+1)
    x = -cos(theta)

    # Compute the C-C weights.
    c = zeros(n+1)
    c[[0,n]] = 1/(n**2-1)
    v = ones(n-1)
    for k in range(1,int(n/2)):
        v -= 2*cos(2*k*theta[1:-1])/(4*k**2-1)
    v -= cos(n*theta[1:-1])/(n**2-1)
    c[1:-1] = 2*v/n

    # Evaluate integrand and integral.
    I = dot(c,f(x))   # use vector inner product
    return I,x

def glint(f,n):
    """
    glint(f,n)

    Perform Gauss-Legendre integration for the function `f` on `n` nodes in (-1,1). Return
    integral and a vector of the nodes used.
    """
    # Nodes and weights are found via a tridiagonal eigenvalue problem.
    beta = 0.5/sqrt(1-(2.0*arange(1,n))**(-2))
    T = diag(beta,1) + diag(beta,-1)
    ev,V = eig(T)
    ev = real_if_close(ev)
    p = argsort(ev)
    x = ev[p]               # nodes
    c = 2*V[0,p]**2         # weights

    # Evaluate the integrand and compute the integral.
    I = dot(c,f(x))         # vector inner product
    return I,x

def intde(f,h,M):
    """
    intde(f,h,M)

    Perform doubly-exponential integration of function `f` over (-Inf,Inf), using
    discretization size `h` and truncation point `M`. Return integral and a vector of the
    nodes used.
    """
    # Find where to truncate the trapezoid sum.
    K = ceil( log(4/pi*log(2*M))/h )

    # Integrate by trapezoids in a transformed variable t.
    t = h*arange(-K,K+1)
    x = sinh(pi/2*sinh(t))
    dxdt = pi/2*cosh(t)*cosh(pi/2*sinh(t))

    I = h*dot(f(x),dxdt)
    return I,x

def intsing(f,h,delta):
    """
    intsing(f,h,delta)

    Integrate function `f` (possibly singular at 1 and -1) over [-1+`delta`,1-`delta`]
    using discretization size `h`. Return integral and a vector of the
    nodes used.
    """
    # Find where to truncate the trapezoid sum.
    K = ceil(log(-2/pi*log(delta/2))/h)

    # Integrate over a transformed variable.
    t = h*arange(-K,K+1)
    x = tanh(pi/2*sinh(t))
    dxdt = pi/2*cosh(t) / (cosh(pi/2*sinh(t))**2)

    I = h*dot(f(x),dxdt)
    return I,x
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from scipy.linalg import *
from numpy.linalg import *
from scipy.optimize import root_scalar
from scipy.integrate import solve_ivp
from FNC04 import levenberg

def shoot(phi,xspan,lval,lder,rval,rder,init):
	"""
	shoot(phi,xspan,lval,lder,rval,rder,init)

	Use the shooting method to solve a two-point boundary value problem. The ODE is
	u'' = `phi`(x,u,u') for x in `xspan`. Specify a function value or derivative at
	the left endpoint using `lval` and `lder`, respectively, and similarly for the
	right endpoint  using `rval` and `rder`. (Use an empty array to denote an
	unknown quantity.) The value `init` is an initial guess for whichever value is
	missing at the left endpoint.

	Return vectors for the nodes, the values of u, and the values of u'.
	"""

	# Tolerances for IVP solver and rootfinder.
	ivp_opt = 1e-6
	optim_opt = 1e-5

	# Evaluate the difference between computed and target values at x=b.
	def objective(s):
		nonlocal x, v   # change these values in outer scope
		# Combine s with the known left endpoint value.
		if len(lder)==0:
			v_init = [ lval[0], s ]  
		else: 
			v_init = [ s, lder[0] ]

		# ODE posed as a first-order equation in 2 variables.
		def shootivp(x,v):
			return array([ v[1], phi(x,v[0],v[1]) ])

		x = linspace(xspan[0],xspan[1],400)  # make decent plots on return
		sol = solve_ivp(shootivp,xspan,v_init,t_eval=x)
		x = sol.t;  v = sol.y

		if len(rder)==0:
			return v[0,-1] - rval[0] 
		else:
			return v[1,-1] - rder[0]

	# Find the unknown quantity at x=a by rootfinding.
	x = [];  v = [];   # the values will be overwritten
	s = root_scalar(objective,x0=init,x1=init+0.05,xtol=optim_opt).root

	# Don't need to solve the IVP again. It was done within the
	# objective function already.
	u = v[0]            # solution
	dudx = v[1]         # derivative

	return x,u,dudx


def diffmat2(n,xspan):
	"""
	diffmat2(n,xspan)

	Compute 2nd-order-accurate differentiation matrices on `n`+1 points in the
	interval `xspan`. Return a vector of nodes, and the matrices for the first
	and second derivatives.
	"""
	a,b = xspan
	h = (b-a)/n
	x = linspace(a,b,n+1)   # nodes

	# Define most of Dx by its diagonals.
	dp = 0.5/h*ones(n)         # superdiagonal
	dm = -0.5/h*ones(n)        # subdiagonal
	Dx = diag(dm,-1) + diag(dp,1)

	# Fix first and last rows.
	Dx[0,:3] = array([-1.5,2,-0.5])/h
	Dx[-1,-3:] = array([0.5,-2,1.5])/h

	# Define most of Dxx by its diagonals.
	d0 =  -2/h**2*ones(n+1)    # main diagonal
	dp =  ones(n)/h**2         # superdiagonal and subdiagonal
	Dxx = diag(d0,0) + diag(dp,-1) + diag(dp,1)

	# Fix first and last rows.
	Dxx[0,:4] = array([2,-5,4,-1])/h**2
	Dxx[-1,-4:] = array([-1,4,-5,2])/h**2

	return x,Dx,Dxx

def diffcheb(n,xspan):
	"""
	diffcheb(n,xspan)

	Compute Chebyshev differentiation matrices on `n`+1 points in the
	interval `xspan`. Return a vector of nodes, and the matrices for the first
	and second derivatives.
	"""
	x = -cos( arange(n+1)*pi/n )   # nodes in [-1,1]
	Dx = zeros([n+1,n+1])
	c = hstack([2.,ones(n-1),2.])    # endpoint factors

	# Off-diagonal entries
	Dx = zeros([n+1,n+1])
	for i in range(n+1):
		for j in range(n+1):
			if i != j:
				Dx[i,j] = (-1)**(i+j) * c[i] / (c[j]*(x[i]-x[j])) 

	# Diagonal entries by the "negative sum trick"
	for i in range(n+1):
		Dx[i,i] = -sum( [Dx[i,j] for j in range(n+1) if j!=i] )    

	# Transplant to [a,b]
	a,b = xspan
	x = a + (b-a)*(x+1)/2
	Dx = 2*Dx/(b-a)

	# Second derivative
	Dxx = Dx @ Dx

	return x,Dx,Dxx

def bvplin(p,q,r,xspan,lval,rval,n):
	"""
		bvplin(p,q,r,xspan,lval,rval,n)

	Use finite differences to solve a linear bopundary value problem. The ODE is
	u''+`p`(x)u'+`q`(x)u = `r`(x) on the interval `xspan`, with endpoint function
	values given as `lval` and `rval`. There will be `n`+1 equally spaced nodes,
	including the endpoints.

	Return vectors of the nodes and the solution values.
	"""
	x,Dx,Dxx = diffmat2(n,xspan)

	P = diag(p(x))
	Q = diag(q(x))
	L = Dxx + P@Dx + Q     # ODE expressed at the nodes

	# Replace first and last rows using boundary conditions.
	I = eye(n+1)
	A = vstack([ I[0], L[1:-1], I[-1] ] )
	b = hstack([ lval, r(x[1:-1]), rval ])

	# Solve the system.
	u = solve(A,b)

	return x,u

def bvp(phi,xspan,lval,lder,rval,rder,init):
	"""
	bvp(phi,xspan,lval,lder,rval,rder,init)

	Use finite differences to solve a two-point boundary value problem. The ODE is
	u'' = `phi`(x,u,u') for x in `xspan`. Specify a function value or derivative at
	the left endpoint using `lval` and `lder`, respectively, and similarly for the
	right endpoint  using `rval` and `rder`. (Use an empty array to denote an
	unknown quantity.) The value `init` is an initial guess for whichever value is
	missing at the left endpoint.

	Return vectors for the nodes and the values of u.
	"""
	n = len(init) - 1
	x,Dx,Dxx = diffmat2(n,xspan)
	h = x[1]-x[0]

	def residual(u):
		# Compute the difference between u'' and phi(x,u,u') at the
		# interior nodes and appends the error at the boundaries.
		dudx = Dx@u                   # discrete u'
		d2udx2 = Dxx@u                # discrete u''
		f = d2udx2 - phi(x,u,dudx)

		# Replace first and last values by boundary conditions.
		if len(lder)==0:
			f[0] = (u[0] - lval[0])/h**2
		else:
			f[0] = (dudx[0] - lder[0])/h
		if len(rder)==0:
			f[-1] = (u[-1] - rval[0])/h**2
		else:
			f[-1] = (dudx[-1] - rder[0])/h
		return f

	u = levenberg(residual,init)
	return x,u[:,-1]

def fem(c,s,f,a,b,n):
	"""
	fem(c,s,f,a,b,n)

	Use a piecewise linear finite element method to solve a two-point boundary
	value problem. The ODE is (`c`(x)u')' + `s`(x)u = `f`(x) on the interval
	[`a`,`b`], and the boundary values are zero. The discretization uses `n` equal
	subintervals.

	Return vectors for the nodes and the values of u.
	"""
	# Define the grid.
	h = (b-a)/n
	x = linspace(a,b,n+1)

	# Templates for the subinterval matrix and vector contributions.
	Ke = array( [[1,-1], [-1,1]] )
	Me = (1/6)*array( [[2,1], [1,2]] )
	fe = (1/2)*array([1, 1])

	# Evaluate coefficent functions and find average values.
	cval = c(x);   cbar = (cval[:-1]+cval[1:]) / 2
	sval = s(x);   sbar = (sval[:-1]+sval[1:]) / 2
	fval = f(x);   fbar = (fval[:-1]+fval[1:]) / 2

	# Assemble global system, one interval at a time.
	K = zeros([n-1,n-1]);  M = zeros([n-1,n-1]);  f = zeros(n-1)
	K[0,0] = cbar[0]/h;  M[0,0] = sbar[0]*h/3;  f[0] = fbar[0]*h/2
	K[-1,-1] = cbar[-1]/h;  M[-1,-1] = sbar[-1]*h/3;  f[-1] = fbar[-1]*h/2
	for k in range(1,n-1):
		K[k-1:k+1,k-1:k+1] += (cbar[k]/h) * Ke
		M[k-1:k+1,k-1:k+1] += (sbar[k]*h) * Me
		f[k-1:k+1] += (fbar[k]*h) * fe

	# Solve system for the interior values.
	u = solve(K+M,f)
	u = hstack([0, u, 0])      # put the boundary values into the result

	return x,u
from scipy import *
from numpy import *
from scipy.linalg import *
from numpy.linalg import *

def diffper(n,xspan):
	"""
	diffper(n,xspan)

	Construct 2nd-order differentiation matrices for functions with periodic end
	conditions, using `n` unique nodes in the interval `xspan`. Return a vector of
	nodes and the  matrices for the first and second derivatives.
	"""
	a,b = xspan
	h = (b-a)/n
	x = a + h*arange(n)   # nodes, omitting the repeated data

	# Construct Dx by diagonals, then correct the corners.
	dp = 0.5/h*ones(n-1)        # superdiagonal
	dm = -0.5/h*ones(n-1)       # subdiagonal
	Dx = diag(dm,-1) + diag(dp,1)
	Dx[0,-1] = -1/(2*h)
	Dx[-1,0] = 1/(2*h)

	# Construct Dxx by diagonals, then correct the corners.
	d0 =  -2/h**2*ones(n)         # main diagonal
	dp =  ones(n-1)/h**2         # superdiagonal and subdiagonal
	Dxx = diag(d0) + diag(dp,-1) + diag(dp,1)
	Dxx[0,-1] = 1/(h**2)
	Dxx[-1,0] = 1/(h**2)

	return x,Dx,Dxx
from scipy import *
from numpy import *
from scipy.linalg import *
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import warnings
from FNC10 import diffmat2

def rectdisc(m,xspan,n,yspan):
	"""
	rectdisc(m,xspan,n,yspan)

	Create matrices and helpers for finite-difference discretization of a rectangle that is
	the tensor  product of intervals `xspan` and `yspan`, using `m`+1 and `n`+1 points in
	the two coordinates.
	"""
	# Initialize grid and finite differences.
	x,Dx,Dxx = diffmat2(m,xspan)
	y,Dy,Dyy = diffmat2(n,yspan)
	X,Y = meshgrid(x,y)

	# Locate boundary points.
	isbndy = tile(True,(n+1,m+1))
	isbndy[1:-1,1:-1] = False

    # Get the diff. matrices recognized as sparse. Also include reshaping functions.
	disc = {
		"Dx":sp.lil_matrix(Dx), "Dxx":sp.lil_matrix(Dxx),
		"Dy":sp.lil_matrix(Dy), "Dyy":sp.lil_matrix(Dyy),
		"Ix":sp.eye(m+1,format="lil"), "Iy":sp.eye(n+1,format="lil"),
		"isbndy":isbndy,
		"vec": lambda U: U.flatten(),
		"unvec": lambda u: reshape(u,(n+1,m+1))
	}
	return X,Y,disc

def poissonfd(f,g,m,xspan,n,yspan):
	"""
	poissonfd(f,g,m,xspan,n,yspan)

	Solve Poisson's equation on a rectangle by finite differences. Function `f` is the
	forcing function and function `g` gives the  Dirichlet boundary condition. The rectangle
	is the tensor product of intervals `xspan` and `yspan`,  and the discretization uses
	`m`+1 and `n`+1 points in the two coordinates.

	Return matrices of the solution values, and the coordinate functions, on the grid.
	"""
	# Initialize the rectangle discretization.
	X,Y,d = rectdisc(m,xspan,n,yspan)

	# Form the collocated PDE as a linear system.
	A = sp.kron(d["Iy"],d["Dxx"]) + sp.kron(d["Dyy"],d["Ix"])  # Laplacian matrix
	b = d["vec"](f(X,Y))

	# Replace collocation equations on the boundary.
	scale = amax(abs(A[n+1,:]))
	I = sp.lil_matrix(sp.eye((m+1)*(n+1)))
	isbndy = d["isbndy"]
	vec = d["vec"]
	A[vec(isbndy),:] = scale*I[vec(isbndy),:]                 # Dirichet assignment
	b[vec(isbndy)] = scale*g( X[isbndy],Y[isbndy] )  # assigned values

	# Solve the linear sytem and reshape the output.
	u = spsolve(A,b)
	U = d["unvec"](u)
	return U,X,Y


def newtonpde(f,g,m,xspan,n,yspan):
	"""
	newtonpde(f,g,m,xspan,n,yspan)

	Newton's method with finite differences to solve the PDE `f`(u,x,y,disc)=0 on the
	rectangle `xspan` ``\times`` `yspan`, subject to `g`(x,y)=0 on the boundary. Use `m`+1
	points in x by `n`+1 points in y.

	Return matrices of the solution values, and the coordinate functions, on the grid.
	"""
	# Discretization.
	X,Y,d = rectdisc(m,xspan,n,yspan)

	# This evaluates the discretized PDE and its Jacobian, with all the
	# boundary condition modifications applied.
	bndy = d["isbndy"]
	vec = d["vec"]
	def residual(U):
		R,J = f(U,X,Y,d)
		scale = amax(abs(J))
		I = sp.eye((m+1)*(n+1),format="lil")
		J[vec(bndy),:] = scale*I[vec(bndy),:]
		XB = X[bndy];  YB = Y[bndy];
		R[bndy] = scale*(U[bndy] - g(XB,YB))
		r = vec(R)
		return r,J

	# Intialize the Newton iteration.
	U = zeros(X.shape)
	r,J = residual(U)
	tol = 1e-10;  itermax = 20;
	s = 2;  normr = norm(r);  k = 1;

	lamb = 1
	I = sp.eye((m+1)*(n+1))
	while (norm(s) > tol) and (normr > tol):
		s = -spsolve(J.T@J + lamb*I, J.T@r)  # damped step
		Unew = U + d["unvec"](s)
		rnew,Jnew = residual(Unew)

		if norm(rnew) < normr:
			# Accept and update.
			lamb = lamb/6;   # dampen the Newton step less
			U = Unew;  r = rnew;  J = Jnew;
			normr = norm(r)
			k = k+1
			print(f"Norm of residual = {normr:.4g}")
		else:
			# Reject.
			lamb = lamb*4;   # dampen the Newton step more

		if k==itermax:
			warnings.warn("Maximum number of Newton iterations reached.")
			break

	return U,X,Y

