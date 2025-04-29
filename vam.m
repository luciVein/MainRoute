cost= input('Enter Cost matrix');

supply = input('enter supply');
demand = input('enter demand');

sum_supp = sum(supply);
sum_dem = sum(demand);

[r,c] = size(cost);

if sum_supp ~= sum_dem

    if sum_supp > sum_dem
        temp = zeros(r,1);
        cost= [cost temp];
        demand= [demand sum_supp-sum_dem];
    else
         temp= zeros(1,c);
         cost = [cost; temp];
         supply = [supply sum_dem-sum_supp] ;
    end

end

[r,c] = size(cost);

disp(supply);
disp(demand); 
disp(cost);

alloc= ones(r,c)*inf;

demo = cost;

% VAM

row_penalty = zeros(1,r);
col_penalty = zeros(1,c);

for i=1:r
    rrow = cost(i,:);
    rrow = sort(rrow);

    if rrow(2) ~= inf
        row_penalty(i) = rrow(2) - rrow(1);
    else
        row_penalty(i) = rrow(1);
    end

end

for j=1:c
    rcol = cost(:,j);
    rcol = sort(rcol);
    
    if rcol(2) ~= inf
        col_penalty(j) = rcol(2) - rcol(1);
    else
        col_penalty(j) = rcol(1);
    end
end


while any(supply>0) || any(demand>0)

    [max_row, row_Idx] = max(row_penalty);
    [max_col, col_Idx] = max(col_penalty);

    % Update allocation/////////////////////////////////
    if max_row > max_col  %  % row penalty is the highest
        [min_cost, col_min_Idx] = min(cost(row_Idx,:));
        
        supply_idx = row_Idx;
        demand_idx = col_min_Idx;

        reduce = min( supply(supply_idx), demand(demand_idx) );
        supply(supply_idx) = supply(supply_idx) - reduce;
        demand(demand_idx) = demand(demand_idx) - reduce;

        alloc(supply_idx, demand_idx) = reduce;

    else
        [min_cost, row_min_Idx] = min(cost(:,col_Idx));
        
        supply_idx = row_min_Idx;
        demand_idx = col_Idx;

        reduce = min( supply(supply_idx), demand(demand_idx) );
        supply(supply_idx) = supply(supply_idx) - reduce;
        demand(demand_idx) = demand(demand_idx) - reduce;

        alloc(supply_idx, demand_idx) = reduce;
    end
    % allocation end

    % fill the lines

    if supply(supply_idx) == 0
        cost(supply_idx,:) = inf;
        row_penalty(supply_idx) = -inf;
        % calculating col penalty again
        fprintf("row is going to be delete\n");

        for j=1:c
            if col_penalty(j) ~= -inf
                rcol = cost(:,demand_idx);
                rcol = sort(rcol);

                if rcol(2) ~= inf
                    col_penalty(j) = rcol(2) - rcol(1);
                else 
                    col_penalty(j) = rcol(1);
                end

            end
        end



    else
        cost(:,demand_idx) = inf;
        col_penalty(demand_idx) = -inf;

        % calculating row penalty again
        fprintf("col is going to be delete\n");

       for i=1:r
            if row_penalty(i) ~= -inf

                rrow = cost(i,:);
                rrow = sort(rrow);

                if rrow(2) ~= inf
                   row_penalty(i) = rrow(2) - rrow(1);
                else
                   row_penalty(i) = rrow(1);
                end
            end
        end

    end


cost

end


ans = 0;
for i = 1 : m
    for j = 1 : n 
        if allocation(i,j) ~= inf
            ans = ans + allocation(i,j) * demo(i,j);
        end 
    end 
end 
allocation
fprintf("ans is %d" , ans );

