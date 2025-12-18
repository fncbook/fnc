function cmap = redsblues()
% Create a red->gray->blue colormap.    
    m = 256;
    m2 = floor(m / 2);
    cmap = zeros(m, 3);

    % Hue
    cmap(1:m2, 1) = 0.95;
    cmap(m2+1:m, 1) = 2/3;
    
    % Saturation
    cmap(1:m2, 2) = linspace(1, 0, m2);
    cmap(m2+1:m, 2) = linspace(0, 1, m-m2);
    
    % Value
    cmap(1:m2, 3) = linspace(0.5, 1, m2);
    cmap(m2+1:m, 3) = linspace(1, 0.5, m-m2);
    
    cmap = hsv2rgb(cmap);

    % Hue
    t = linspace(-1, 1, m)';
    cmap(:, 1) = 240*(tanh(4*t)+1)/720;
    
    % Saturation
    cmap(1:m2, 2) = linspace(1, 0, m2);
    cmap(m2+1:m, 2) = linspace(0, 1, m-m2);
    
    % Value
    cmap(1:m2, 3) = linspace(0.4, 1, m2);
    cmap(m2+1:m, 3) = linspace(1, 0.4, m-m2);
    
    cmap = hsv2rgb(cmap);

end