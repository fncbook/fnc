function du_dt = predprey(t, u, p)
    alpha = p(1);  beta = p(2);
    y = u(1);      z = u(2);
    s = (y * z) / (1 + beta * y);  % appears in both equations
    du_dt = [ y * (1 - alpha * y) - s;  -z + s ];
end
