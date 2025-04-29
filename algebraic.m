function x = algebraic(A, b, cT)
    [m, n] = size(A);
    combinations = nchoosek(1:n, m);
    z = -Inf;

    for y = 1: size(combinations, 1)
        B = zeros(m, m);
        for i = 1: m
            B(:, i) = A(:, combinations(y, i));
        end
        
        xt = zeros(n, 1);
        if det(B) == 0
            continue;
        end
        xb = B\b; % xb = inverse(B) * b

        for i = 1: m
            xt(combinations(y, i)) = xb(i);
        end
        argval = cT*xt;
        if argval > z
            possible = 1;
            for i = 1: n
                if xt(i) < 0
                    possible = 0;
                    break
                end
            end

            if possible == 1
                z = argval
                x = xt;
            end
        end
    end
end
