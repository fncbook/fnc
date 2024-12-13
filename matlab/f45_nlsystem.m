function [f, J] = nlsystem(x)
    f = zeros(3, 1);   % ensure a column vector output
    f(1) = exp(x(2) - x(1)) - 2;
    f(2) = x(1) * x(2) + x(3);
    f(3) = x(2) * x(3) + x(1)^2 - x(2);
    J(1, :) = [-exp(x(2) - x(1)), exp(x(2) - x(1)), 0];
    J(2, :) = [x(2), x(1), 1];
    J(3, :) = [2 * x(1), x(3)-1, x(2)];
end