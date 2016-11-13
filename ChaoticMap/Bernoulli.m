function [ randB ] = Bernoulli(rows, cols)
% where a is a control parameter in (0, 1)
a = 0.5;
x0 = @rand;
x0 = feval(x0);

switch nargin
    case 1
        randB = zeros(rows, rows);
        for j = 1 : rows
            for i = 1 : rows
                if 0 <= x0.x0 <= a
                    randB(i,j) = x0.x0 / (1 - a);
                end
                if 1-a <= x0.x0 <= 1
                    randB(i,j) = (x0.x0 - (1 - a)) / a;
                end
                x0.x0 = feval(x0);
            end
        end
    case 2
        randB = zeros(rows, cols);
        for j = 1 : cols
            for i = 1 : rows
                if 0 <= x0.x0 <= a
                    randB(i,j) = x0.x0 / (1 - a);
                end
                if 1-a <= x0.x0 <= 1
                    randB(i,j) = (x0.x0 - (1 - a)) / a;
                end
                x0.x0 = feval(x0);
            end
        end
    otherwise
        if 0 <= x0.x0 <= a
            randB = x0.x0 / (1 - a);
        end
        if 1-a <= x0.x0 <= 1
            randB = (x0.x0 - (1 - a)) / a;
        end
end

end

