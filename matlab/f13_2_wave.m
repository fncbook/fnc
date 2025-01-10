function dw_dt = timederiv(t, w, p)
    [Dxx, Dyy, pack, unpack] = p{:};
    [U, V] = unpack(w);
    dU_dt = V;
    dV_dt = Dxx * U + U * Dyy';
    dw_dt = pack(dU_dt, dV_dt);
end