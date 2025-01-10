function udot = pendulums(t, u, p)
    gamma = p(1);  L = p(2);  k = p(3);
    g = 9.8;
    udot = zeros(4, 1);
    udot(1:2) = u(3:4);
    udot(3) = -gamma * u(3) - (g / L) * sin(u(1)) + k * (u(2) - u(1));
    udot(4) = -gamma * u(4) - (g / L) * sin(u(2)) + k * (u(1) - u(2));
end
