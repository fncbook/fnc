# Glossary

:::{glossary}
adjacency matrix
: Matrix whose nonzero entries show the links between nodes in a graph. (@definition-adjacency-matrix)

adjoint
: Conjugate transpose of a complex matrix. (@definition-adjoint)

advection equation
: Archetypical PDE of hyperbolic type, representing transport phenomena. (@definition-advection-equation)

algorithm
: List of instructions for transforming data into a result. (@section-intro-algorithms)

Arnoldi iteration
: Stable algorithm for finding orthonormal bases of nested Krylov subspaces. (@definition-arnoldi-iteration)

asymptotic
: Relationship indicating that two functions have the same leading behavior in some limit. (@definition-asymptotic-notation)

backward error
: Change to the input of a problem required to produce the result found by an inexact algorithm. (@definition-backward-error)

backward substitution
: Systematic method for solving a linear system with an upper triangular matrix. @eq-backsub

bandwidth
: The number of diagonals around the main diagonal that have nonzero elements. (@definition-bandwidth)

barycentric interpolation formula
: Computationally useful expression for the interpolating polynomial as a ratio of rational terms. @bary2

big-O
: Relationship indicating that one function is bounded above by a multiple of another in some limit. (@definition-asymptotic-notation)

boundary-value problem
: A differential equation with which partial information about the solution is given at multiple points on the boundary of the domain. (@definition-tpbvp)

cardinal function
: Interpolating function that is 1 at one node and 0 at all the others. (@definition-cardinal-function)

Cholesky factorization
: Symmetrized version of LU factorization for SPD matrices. (@theorem-cholesky-factorization)

collocation
: Solution of a differential equation by imposing it approximately at a set of nodes. @fdlinnobc

condition number
: Ratio of the size of change in the output of a function to the size of change in the input that produced it. (@definition-condition-number)

cubic spline
: Piecewise cubic function with two globally continuous derivatives, most often used for interpolation or approximation. (@definition-cubic-spline)

diagonalizable matrix
: Matrix that admits an eigenvalue decomposition. Also known as *nondefective*. (@definition-evd)

differentiation matrix
: Matrix mapping a vector of function values to a vector of approximate derivative values. (@definition-diffmat)

Dirichlet condition
: Boundary condition specifying the value of the solution. (@definition-bctypes)

dominant eigenvalue
: Eigenvalue with the largest magnitude (absolute value, in the real case). (@definition-dominant-eigenvalue)

double precision
: Typical standard in floating-point representation, using 64 bits to achieve about 16 decimal significant digits of precision. @doubleprec

eigenvalue
: Scalar $\lambda$ such that $\mathbf{A}\mathbf{x} = \lambda \mathbf{x}$ for a square matrix $\mathbf{A}$ and nonzero vector $\mathbf{x}$. (@definition-eigenvalue)

eigenvalue decomposition
: Expression of a square matrix as the product of eigenvector and diagonal eigenvalue matrices. Abbreviated EVD. (@definition-evd)

eigenvector
: Vector for which the action of a matrix is effectively one-dimensional. (@definition-eigenvalue)

Euler's method
: Prototype of all IVP solution methods, obtained by assuming constant derivatives for the solution over short time intervals. (@definition-eulerivp)

evolutionary PDE
: A partial differential equation in which one of the independent variables is time or a close analog. (@definition-evolutionary-pde)

extrapolation
: Use of multiple discretization values to cancel out leading terms in an error expansion. (@section-integration-extrapolation)

finite-difference formula
: Linear combination of function values that approximates the value of a derivative of the function at a point. (@definition-finite-difference)

finite element method (FEM)
: Use of piecewise integration to pose a linear system of equations for the approximate solution of a boundary-value problem. (@section-galerkin-fem)

fixed-point iteration
: Repeated application of a function in hopes of converging to a fixed point. (@definition-fpiteration)

fixed-point problem
: Finding a value of a given function where the input and output values are the same; equivalent to rootfinding. (@definition-fpproblem)

floating-point numbers
: A finite set that substitutes for the real numbers in machine calculations. Denoted by $\mathbb{F}$. (@definition-float)

flop
: A single arithmetic operation on floating-point numbers, often counted as a proxy for computer runtime. (@definition-flop)

forward substitution
: Systematic method for solving a linear system with a lower triangular matrix. @forwardsub

Frobenius norm
: Matrix norm computed by applying the vector 2-norm to a vector interpretation of the matrix. (@definition-frobenius)

Gauss–Newton method
: Generalization of Newton's method for nonlinear least squares. (@section-nonlineqn-nlsq)

Gaussian elimination
: Use of row operations to transform a linear system to an equivalent one in triangular form. (@section-linsys-lu)

generating polynomials
: A pair of polynomials whose coefficients match those of a multistep method for IVPs. (@section-ivp-multistep)

global error
: Error made by an IVP method over the entire time interval of the solution. (@section-ivp-euler)

GMRES
: Iterative solution of a linear system through stable least-squares solutions on nested Krylov subspaces. (@section-krylov-gmres)

graph
: Representation of a network as a set of nodes and edges. (@section-matrixanaly-insight)

hat functions
: Cardinal functions for piecewise linear interpolation. (@section-localapprox-pwlin)

heat equation
: Archetypical parabolic PDE that describes diffusion. (@section-diffusion-blackscholes)

hermitian
: Combination of transpose and elementwise complex conjugation. Same as the {term}`adjoint`.(@definition-adjoint)

hermitian matrix
: A matrix that equals its own adjoint. Also known as _self-adjoint_.(@definition-symmetric-matrix)

hermitian positive definite matrix
: Matrix that is hermitian with strictly positive eigenvalues; complex variant of symmetric positive definite. (Abbreviated HPD matrix.)(@section-matrixanaly-symm-eig)

homogeneous boundary condition
: Having a zero value. (@section-bvp-tpbvp, @section-diffusion-blackscholes)

identity matrix
: Matrix with ones on the diagonal and zeros elsewhere, acting as the multiplicative identity. (@definition-identity-matrix)

ill-conditioned
: Exhibiting a large condition number, indicating high sensitivity of a result to changes in the data. (@section-intro-conditioning)

implicit
: Formula that defines a quantity only indirectly, e.g., as the solution of a nonlinear equation. (@section-ivp-multistep)

induced matrix norm
: Norm computed using the interpretation of a matrix as a linear operator. (@section-linsys-norms)

initial-value problem
: An ordinary differential equation (possibly vector-valued) together with an initial condition. (@section-ivp-basics, @section-ivp-systems)

inner product
: Scalar or dot product of a pair of vectors, or its extension to a pair of functions. (@section-globalapprox-orthogonal)

interpolation
: Construction of a function that passes through a given set of data points. (@section-linsys-polyinterp, @section-localapprox-interpolation)

inverse iteration
: Subtraction of a shift followed by matrix inversion, used in power iteration to transform the eigenvalue closest to a target value into a dominant one.  (@section-krylov-inviter)

Jacobian matrix
: Matrix of first partial derivatives that defines the linearization of a vector-valued function. (@section-nonlineqn-newtonsys)

Kronecker product
: Alternative type of matrix multiplication useful for problems on a tensor-product domain. (@section-twodim-laplace)

Krylov subspace
: Vector space generated by powers of a square matrix that is often useful for reducing the dimension of large problems. (@section-krylov-subspace)

Lagrange formula
: Theoretically useful expression for an interpolating polynomial. (@section-globalapprox-polynomial)

Lanczos iteration
: Specialization of the Arnoldi iteration to the case of a hermitian (or real symmetric) matrix. (@section-krylov-minrescg)

Laplace equation
: Archetypical elliptic PDE describing a steady state. (@section-twodim-laplace)

linear convergence
: Sequence in which the difference between sequence value and limit asymptotically decreases by a constant factor at each term, making a straight line on a log-linear graph.\ 
(@section-nonlineqn-fixed-point)

linear least-squares problem
: Minimization of the 2-norm of the residual for an overdetermined linear system. (@section-leastsq-fitting)

local truncation error
: Discretization error made in one time step of an IVP solution method. (@section-ivp-euler, @section-ivp-multistep)

LU factorization
: Factorization of a square matrix into the product of a unit lower triangular matrix and an upper triangular matrix. (@section-linsys-lu)

machine epsilon
: Distance from 1 to the next-largest floating-point number. Also called unit roundoff or machine precision, though the usages are not consistent across different references.\ 
(@section-intro-floating-point)

matrix condition number
: Norm of the matrix times the norm of its inverse, equivalent to the condition number for solving a linear system. (@section-linsys-condition-number)

method of lines
: Solution technique for partial differential equations in which each independent variable is discretized separately. (@section-diffusion-methodlines)

multistep
: Formula using information over more than a single time step to advance the solution. (@section-ivp-multistep)

Neumann condition
: Boundary condition specifying the derivative of the solution. (@section-bvp-tpbvp)

Newton's method
: Rootfinding iteration that uses the linearization of the given function in order to define the next root approximation. (@section-nonlineqn-newton)

nodes
: Values of the independent variable where an interpolant's values are prescribed. (@section-localapprox-interpolation)

nonlinear least-squares problem
: Minimization of the 2-norm of the residual of a function that depends nonlinearly on a vector. (@section-nonlineqn-nlsq)

norm
: Means of defining the magnitude of a vector or matrix. (@section-linsys-norms)

normal matrix
: Matrix that has a {term}`unitary matrix` of eigenvectors in an {term}`eigenvalue decomposition`. (@section-matrixanaly-evd)

normal equations
: Square linear system equivalent to the linear least-squares problem. (@section-leastsq-normaleqns)

numerical integration
: Estimation of a definite integral by combining values of the integrand, rather than by finding an antiderivative. (@section-localapprox-integration)

one-step IVP method
: IVP solver that uses information from just one time level to advance to the next. (@section-ivp-euler)

ONC matrix
: Matrix whose columns are orthonormal vectors. (@section-leastsq-qr)

order of accuracy
: Leading power of the truncation error as a function of a discretization size parameter. (@section-localapprox-fd-converge, @section-localapprox-integration, @section-ivp-euler, @section-ivp-multistep)

orthogonal vectors
: Nonzero vectors that have an inner product of zero. (@section-leastsq-qr)

orthogonal matrix
: Square ONC matrix, i.e., matrix whose transpose is its inverse. (@section-leastsq-qr)

orthogonal polynomials
: Family of polynomials whose distinct members have an integral inner product equal to zero, as with Legendre and Chebyshev polynomials. (@section-globalapprox-orthogonal)

orthonormal vectors
: Vectors that are both mutually orthogonal and all of unit 2-norm. (@section-leastsq-qr)

outer product
: Multiplication of two vectors resulting in a rank-1 matrix. (@section-linsys-lu)

overdetermined
: Characterized by having more constraints than available degrees of freedom. (@section-leastsq-fitting)

piecewise linear
: Function that is linear between each consecutive pair of nodes, but whose slope may jump at the nodes. (@section-localapprox-pwlin)

PLU factorization
: LU factorization with row pivoting. (@section-linsys-pivoting)

power iteration
: Repeated application of a matrix to a vector, followed by normalization, resulting in convergence to an eigenvector for the dominant eigenvalue. (@section-krylov-power)

preconditioning
: Use of an approximate inverse to improve the convergence rate of Krylov iterations for a linear system. (@section-krylov-precond)

pseudoinverse
: Rectangular matrix that maps data to solution in the linear least-squares problem, generalizing the matrix inverse. (@section-leastsq-normaleqns)

QR factorization
: Representation of a matrix as the product of an orthogonal and an upper triangular matrix. (@section-leastsq-qr)

quadratic convergence
: Sequence in which the difference between sequence value and limit asymptotically decreases by a constant times the square of the preceding difference. (@section-nonlineqn-newton)

quasi-Newton methods
: Rootfinding methods that overcome the issues of Jacobian computation and lack of global convergence in Newton's method. (@section-nonlineqn-quasinewton)

quasimatrix
: Collection of functions (such as orthogonal polynomials) that have algebraic parallels to columns of a matrix. (@section-globalapprox-orthogonal)

Rayleigh quotient
: Function of vectors that equals an eigenvalue when given its eigenvector as input. (@section-matrixanaly-symm-eig)

reduced QR factorization
: See *thin QR*.

reduced SVD
: See *thin SVD*.

residual
: For a linear system, the difference between $\mathbf{b}$ and $\mathbf{A}\tilde{\mathbf{x}}$ for a computed solution approximation $\tilde{\mathbf{x}}$. More generally, the actual value of a quantity that is made zero by an exact solution. (@section-linsys-condition-number, @section-nonlineqn-rootproblem)

restarting
: Technique used in GMRES to prevent the work per iteration and overall storage from growing uncontrollably. (@section-krylov-gmres)

rootfinding problem
: Finding the input value for a given function which makes that function zero. (@section-nonlineqn-rootproblem)

row pivoting
: Reordering rows during LU factorization to ensure that the factorization exists and can be computed stably. (@section-linsys-pivoting)

Runge phenomenon
: Manifestation of the instability of polynomial interpolation at equally spaced nodes as degree increases. (@section-globalapprox-stability)

Runge–Kutta
: One-step method for IVPs that evaluates the derivative of the solution more than once to advance a single step. (@section-ivp-rk)

secant method
: Scalar quasi-Newton method that uses a secant line rather than a tangent line to define a root estimate. (@section-nonlineqn-secant)

shooting
: Unstable technique for solving a boundary-value problem in which an initial value is sought for by a rootfinding algorithm. (@section-bvp-shooting)

similar matrices
: Matrices that are linked by a {term}`similarity transformation`, thus sharing the same set of eigenvalues. (@section-matrixanaly-evd)

similarity transformation
: Mapping $\mathbf{A} \mapsto \mathbf{S}\mathbf{A}\mathbf{S}^{-1}$ for an invertible $\mathbf{S}$. The transformation leaves eigenvalues, but not eigenvectors, unchanged. (@section-matrixanaly-evd)

simple root
: Root of a function at which the derivative of the function is nonzero. (@section-nonlineqn-rootproblem)

singular value decomposition
: Expression of a matrix as a product of two orthogonal/unitary matrices and a nonnegative diagonal matrix. (Abbreviated SVD.) (@section-matrixanaly-svd)

sparse matrix
: Describing a matrix that has mostly zero elements for structural reasons. (@section-linsys-structure, @section-krylov-structure)

spectral convergence
: Exponentially rapid decrease in error as the number of interpolation nodes increases, e.g., as observed in Chebyshev polynomial and trigonometric interpolation. (@section-globalapprox-stability)

stability region
: Region of the complex plane describing when numerical solution of a linear IVP is bounded as $t\to\infty$. (@section-diffusion-absstab)

step size
: Increment in time between successive solution values in a numerical IVP solver. (@section-ivp-euler)

stiff differential equation
: Describes an IVP in which stability is a greater restriction than accuracy for many solution methods, usually favoring the use of an implicit time stepping method. (@section-ivp-implicit, @section-diffusion-stiffness)

subtractive cancellation
: Growth in relative error that occurs when two numbers are added/subtracted to get a result that is much smaller in magnitude than the operands; also called *loss of significance* or *cancellation error*. (@section-intro-conditioning)

superlinear convergence
: Sequence for which the convergence is asymptotically faster than any linear rate. (@definition-superlinear-convergence)

symmetric matrix
: Square matrix that is equal to its transpose. (@definition-symmetric-matrix)

symmetric positive definite matrix
: Matrix that is symmetric and positive definite, thereby permitting a Cholesky factorization. Correspondingly called hermitian positive definite in the complex case. (@section-linsys-structure)

tensor-product domain
: A domain that can be parameterized using variables that lay in a logical rectangle or cuboid; i.e., each variable independently varies in an interval. (@section-twodim-tensorprod)

thin QR factorization
: Variant of the QR factorization that discards information not needed to fully represent the original matrix. (@section-leastsq-qr)

thin SVD
: Variant of the singular value decomposition that discards information not needed to fully represent the original matrix. (@section-matrixanaly-svd)

trapezoid formula
: Numerical integration method resulting from integration of a piecewise linear interpolant. (@section-localapprox-integration, @section-ivp-multistep)

triangular matrix
: Matrix that is all zero either above (for lower triangular) or below (for upper triangular) the main diagonal. (@section-linsys-linear-systems)

tridiagonal matrix
: Matrix with nonzeros only on the main diagonal and the adjacent two diagonals. (@section-linsys-structure)

trigonometric interpolation
: Interpolation of a periodic function by a linear combination of real or complex trigonometric functions. (@section-globalapprox-trig)

truncation error
: Difference between an exact value and an approximation, such as one that truncates an infinite series. (@section-localapprox-fd-converge, @section-localapprox-integration)

unit triangular matrix
: Triangular matrix that has a 1 in every position on the main diagonal. (@section-linsys-lu)

unit vector
: A vector whose norm equals 1. (@section-linsys-norms)

unitary matrix
: Square matrix with complex-valued entries whose columns are orthonormal. (@section-matrixanaly-evd)

unstable
: Allowing perturbations of the data to have much larger effects on the results than can be explained by the problem's condition number. (@section-intro-stability)

upper Hessenberg matrix
: Matrix that has nonzeros only in the upper triangle and first subdiagonal. (@section-krylov-subspace)

Vandermonde matrix
: Matrix whose columns are formed from elementwise powers of a given vector, important for polynomial interpolation and approximation of data. (@section-linsys-polyinterp)

weights
: Coefficients in a linear combination of function values in a finite-difference or integration method. (@section-localapprox-finitediffs, @section-localapprox-integration, @section-globalapprox-integration)

zero-stability
: Boundedness property of multistep methods that is required for convergence. (@section-ivp-zerostability)

:::
