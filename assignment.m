function x = assignment (c)
    
    % Disrespect the 5, I put your - in the food chain.

    n = size(c, 1);
    assert(n==size(c, 2));
    x = zeros(n);
    
    cols = 1:n;
    rows = 1:n;
    % Eliminate smallest element in each row of c
    for j = rows
        c(j, :) = c(j, :) - min(c(j, :));
    end

    % Eliminate smallest element in each coloumn of c
    for i = cols
        c(:, i) = c(: , i) - min(c(:, i));
    end

    crossed = zeros(n);

    round_enrectangled = 1;
    total_enrectangled = 0;

    % Finding largest set of Independent Zeros
    while round_enrectangled ~= 0
        round_enrectangled = 0;
        for j = rows
            num_zeros = 0;
            zero_idx = 0;
            for i = cols
                if c(j, i) == 0 && crossed(j, i) == 0 && x(j, i) == 0
                    num_zeros = num_zeros + 1;
                    zero_idx = i;
                end
            end
            if num_zeros == 1
                round_enrectangled = round_enrectangled + 1;
                x(j, zero_idx) = 1;
                for j2 = rows
                    if c(j2, zero_idx) == 0 && x(j2, zero_idx) == 0
                        crossed(j2, zero_idx) = 1;
                    end
                end
            end
        end

        for i = cols
            num_zeros = 0;
            zero_jdx = 0;
            for j = rows
                if c(j, i) == 0 && crossed(j, i) == 0 && x(j, i) == 0
                    num_zeros = num_zeros + 1;
                    zero_jdx = j;
                end
            end
            if num_zeros == 1
                round_enrectangled = round_enrectangled + 1;
                x(zero_jdx, i) = 1;
                for i2 = cols
                    if c(zero_jdx, i2) == 0 && x(zero_jdx, i2) == 0
                        crossed(zero_jdx, i2) = 1;
                    end
                end
            end
        end
        total_enrectangled = total_enrectangled + round_enrectangled;
    end
    
    if total_enrectangled == n
        return;
    end

    round_marked = 1;
    marked_rows = zeros(n, 1);
    marked_cols = zeros(n, 1);
    
    % Marking all rows for which assignments have not been made.
    for j = rows
        if sum(x(j, :)) == 0 && marked_rows(j) == 0
            marked_rows(j) = 1;
        end
    end
    
    while round_marked ~= 0
        % Marking all columns, not already marked, which have unaasigned zeros in the marked rows.
        round_marked = 0;
        for i = cols
            if marked_cols(i) == 1
                continue
            end
            for j = rows
                if marked_rows(j) == 0
                    continue;
                end
                if crossed(j, i) == 1
                    marked_cols(i) = 1;
                    round_marked = round_marked + 1;
                end
            end
        end

        % Marking all rows, not already marked, which have assigned zeros in the marked coloumns.
        for j = rows
            if marked_rows(j) == 1
                continue
            end
            for i = cols
                if marked_cols(i) == 0
                    continue
                end
                if x(j, i) == 1
                    marked_rows(j) = 1;
                    round_marked = round_marked + 1;
                end
            end
        end

        % Repeating this until no further marking is made.
    end

    col_line = find(marked_cols);
    unmarked_row = find(marked_rows);
    unmarked_col = setdiff(cols, col_line);
    row_line = setdiff(rows, unmarked_row);

    value = min(min(c(unmarked_row, unmarked_col)));
    c(unmarked_row, unmarked_col) = c(unmarked_row, unmarked_col) - value;
    c(row_line, col_line) = c(row_line, col_line) + value;

    % Recursively calling assignment on updated c.
    x = assignment(c);
end
