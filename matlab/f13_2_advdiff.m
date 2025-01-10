function du_dt = timederiv(t, u, p)
    [ep, Dx, Dxx, Dy, Dyy, pack, unpack] = p{:};
    U = unpack(u);
    Uxx = Dxx * U;  Uyy = U * Dyy'; 
    dU_dt = 1 - Dx * U + ep * (Uxx + Uyy);  % PDE
    du_dt = pack(dU_dt);
end