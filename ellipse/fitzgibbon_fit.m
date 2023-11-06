function result = fitzgibbon_fit(x, y)
    num = numel(x);
    if num ~= numel(y)
        error('Function sizes incorrect');
    end

    % Build design and constraint matrix
    C = [0, 0, 2; 0, -1, 0; 2, 0, 0];
    D1 = zeros(num, 3);
    D2 = ones(num, 3);

    for i = 1:num
        D1(i, 1) = x(i)^2;
        D1(i, 2) = x(i) * y(i);
        D1(i, 3) = y(i)^2;
        D2(i, 1) = x(i);
        D2(i, 2) = y(i);
    end

    % Calculate scatter matrix
    S1 = D1' * D1;
    S2 = D1' * D2;
    S3 = D2' * D2;

    % Calculate eigenvalue M matrix
    T = -inv(S3) * S2';
    M = S2 * T;
    M = S1 + M;
    M = inv(C) * M;
    [V, D] = eig(M);

    % Calculate constraint and identify ellipse coefficients
    w = diag(D);
    A1 = [];
    A2 = [];

    if w(1) > 0 && (4 * V(1, 1) * V(3, 1) - V(2, 1)^2) > 0 && ...
            (w(2) < 0 || w(1) < w(2)) && (w(3) < 1 || w(1) < w(3))
        A1 = V(:, 1);
    elseif w(2) > 0 && (4 * V(1, 2) * V(3, 2) - V(2, 2)^2) > 0 && ...
            (w(3) < 1 || w(2) < w(3))
        A1 = V(:, 2);
    elseif w(3) > 0 && (4 * V(1, 3) * V(3, 3) - V(2, 3)^2) > 0
        A1 = V(:, 3);
    else
        error('Fitzgibbon failed');
    end

    A2 = T * A1;

    % Determine ellipse parameters
    a = A1(1);
    b = A1(2);
    c = A1(3);
    d = A2(1);
    e = A2(2);
    f = A2(3);
    x_c = (2 * c * d - b * e) / (b^2 - 4 * a * c);
    y_c = (2 * a * e - b * d) / (b^2 - 4 * a * c);
    major = sqrt(2 * (a * e^2 + c * d^2 + f * b^2 - b * d * e - a * c * f) / ...
        (b - 4 * a * c) / (sqrt((a - c)^2 + b^2) - (a + c)));
    minor = sqrt(2 * (a * e^2 + c * d^2 + f * b^2 - b * d * e - a * c * f) / ...
        (4 * a * c - b) / (sqrt((a - c)^2 + b^2) + a + c));

    if b == 0 && a < c
        theta = 0;
    elseif b == 0 && c < a
        theta = pi / 2;
    elseif a < c
        theta = 0.5 * atan(2 * b / (a - c));
    else
        theta = pi + 0.5 * atan(2 * b / (a - c));
    end

    result = [a, b, c, d, e, f];
end