"""
forwardsub(L,b)

Solve the lower-triangular linear system with matrix `L` and
right-hand side vector `b`.
"""
function forwardsub(L,b)

n = size(L,1)
x = zeros(n)
x[1] = b[1]/L[1,1]
for i = 2:n
    s = sum( L[i,j]*x[j] for j=1:i-1 )
    x[i] = ( b[i] - s ) / L[i,i]
end

return x
end

"""
backsub(U,b)

Solve the upper-triangular linear system with matrix `U` and
right-hand side vector `b`.
"""
function backsub(U,b)

n = size(U,1)
x = zeros(n)
x[n] = b[n]/U[n,n]
for i = n-1:-1:1
    s = sum( U[i,j]*x[j] for j=i+1:n )
    x[i] = ( b[i] - s ) / U[i,i]
end

return x
end

"""
lufact(A)

Compute the LU factorization of square matrix `A`, returning the
factors.
"""
function lufact(A)

n = size(A,1)
L = Matrix(Diagonal(ones(n)))
U = float(copy(A))

# Gaussian elimination
for j = 1:n-1
    for i = j+1:n
        L[i,j] = U[i,j] / U[j,j]   # row multiplier
        U[i,j:n] -= L[i,j]*U[j,j:n]
    end
end

return L,triu(U)
end


function outerlu(A)
  
n = size(A,1); 
Aj = copy(A)
L = zeros(n,n)
U = zeros(n,n)

for j = 1:n
    U[j,j:n] = Aj[j,j:n]
    L[j:n,j] = Aj[j:n,j] / U[j,j]
    Aj[j+1:n,j+1:n] -= L[j+1:n,j]*U[j,j+1:n]'
end

return L,U
end

function lupivot(A)
  
  n = size(A,1); 
  Aj = copy(A)
  L = zeros(n,n)
  U = zeros(n,n)
  pivot = zeros(Int,n)
  unused = trues(n)
  
  for j = 1:n
      pivot[j] = argmax( abs.(Aj[:,j]) )
      U[j,j:n] = Aj[pivot[j],j:n]
      L[unused,j] = Aj[unused,j] / U[j,j]
      unused[pivot[j]] = false
      Aj[unused,j+1:n] -= L[unused,j]*U[j,j+1:n]'
      Aj[pivot[j],j+1:n] .= 0
  end
  
  return L,U,pivot
  end

  function lupivoteasy(A)
  
    n = size(A,1); 
    Aj = copy(A)
    L = zeros(n,n)
    U = zeros(n,n)
    pivot = zeros(Int,n)
    
    for j = 1:n
        pivot[j] = argmax( abs.(Aj[:,j]) )
        U[j,:] = Aj[pivot[j],:]
        L[:,j] = Aj[:,j] / U[j,j]
        Aj -= L[:,j]*U[j,:]'
    end
    
    return tril(L),triu(U),pivot
    end