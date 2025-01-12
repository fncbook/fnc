function dw_dt = wave(t, w, param)
    [c, m, Dx, chop, extend] = param{:};
    u = extend(w(1:m-1));
    z = w(m:2*m);
    du_dt = Dx * z;
    dz_dt = c.^2 .* (Dx * u);
    dw_dt = [ chop(du_dt); dz_dt ];
end