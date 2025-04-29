function x =  dual_simplex(A_leq, A_eq, A_geq, b_leq, b_eq, b_geq, cT)
    % count = size(A, 1) for {leq, eq, geq}
    count_leq = size(A_leq, 1);
    count_eq = size(A_eq, 1);
    count_geq = size(A_geq, 1);

    A_leq = [A_leq eye(count_leq) zeros(count_leq, count_geq+2*count_eq)];
    A_geq = [-A_geq zeros(count_geq, count_leq) eye(count_geq) zeros(count_geq, 2*count_eq)];
    A_eq_pos = [A_eq zeros(count_eq, count_leq+count_geq) eye(count_eq) zeros(count_eq, count_eq)];
    A_eq_neg = [-A_eq zeros(count_eq, count_leq+count_geq+count_eq) eye(count_eq)];
    
    % sz = number of decision variables
    sz = size(cT, 2); 

    A = [A_leq; A_geq; A_eq_pos; A_eq_neg];
    b = [b_leq; -b_geq; b_eq; -b_eq];

    cT = [cT zeros(1, count_leq + count_geq + 2*count_eq)];
    
    [m, n] = size(A);

    v = n-m+1:n;
    xb = b;
    y = A;

    x = zeros(sz, 1);
    
    while true
        y
        xb
        v
        leaving_idx = 0;
        min_xb = 0;
        for j = 1: m
            if xb(j) < min_xb
                min_xb = xb(j);
                leaving_idx = j;
            end
        end

        if leaving_idx == 0
            disp("Optimal solution found.");
            break;
        end

        entering_idx = 0;
        max_ratio = -Inf;
        for i = 1: n
            diff = cT(v) * y(:, i) - cT(i);
            ratio = diff / y(leaving_idx, i);
            if ratio > max_ratio && ratio < 0
                max_ratio = ratio;
                entering_idx = i;
            end
        end

        if entering_idx == 0
            disp("The system has an unbounded solution or no solutions, returning x=0");
            return;
        end
        leaving_idx
        entering_idx

        pivot_value = y(leaving_idx, entering_idx);
        v(leaving_idx) = entering_idx;
        xb(leaving_idx) = xb(leaving_idx) / pivot_value;
        y(leaving_idx, :) = y(leaving_idx, :) / pivot_value;

        for j = 1:m
            if j == leaving_idx
                continue
            end
            factor = y(j, entering_idx) / y(leaving_idx, entering_idx);
            y(j, :) = y(j, :) - y(leaving_idx, :) * factor;
            xb(j) = xb(j) - xb(leaving_idx) * factor;
        end
    end

    for i = 1:m
        if v(i) <= sz
            x(v(i)) = xb(i);
        end
    end

end