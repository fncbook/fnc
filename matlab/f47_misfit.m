function [f, J] = misfit(c, s, w)
    V = c(1);   Km = c(2);
    f = V * s ./ (Km + s) - w;
    J(:,1) = s ./ (Km + s);            % d/d(V)
    J(:,2) = -V * s ./ (Km + s).^2;    % d/d(Km)
end
