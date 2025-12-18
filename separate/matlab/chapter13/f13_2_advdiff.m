function dw_dt = timederiv(t, w, p)
    [ep, Dx, Dxx, Dy, Dyy, pack, unpack] = p{:};
    U = unpack(w);
    Uxx = Dxx * U;  Uyy = U * Dyy'; 
    dU_dt = 1 - Dx * U + ep * (Uxx + Uyy);  % PDE
    dw_dt = pack(dU_dt);
end