function [U, V] = unpack(w, N, unvec_u, unvec_v, extend)
    U = extend( unvec_u(w(1:N)) );
    V = unvec_v(w(N+1:end));
end