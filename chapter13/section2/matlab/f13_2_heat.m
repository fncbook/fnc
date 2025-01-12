function du_dt = timederiv(t, u, p)
    [alpha, Dxx, Dyy, vec, unvec] = p{:};
    U = unvec(u);
    Uxx = Dxx * U;  Uyy = U * Dyy';     % 2nd partials
    dU_dt = alpha * (Uxx + Uyy);  % PDE
    du_dt = vec(dU_dt);
end