function [ Rand ] = Tent(rows, cols, OPTIONS)
OPTIONS_.x0 = @randn;
if  nargin == 3 && ~isempty(OPTIONS) && isstruct(OPTIONS)
    if isfield(OPTIONS, 'x0'); OPTIONS_.x0 = OPTIONS.x0; end
end
x0 = feval(OPTIONS_.x0);

%% Tent function
    function [y] = fcn (x0)
        if x0 < 0.5
            y = 2 * x0;
        end
        if x0 >= 0.5
            y = 2 * (1 - x0);
        end
    end

switch nargin
    case 1
        cols = rows;
        Rand = RCs(rows, cols);
    case 2
        Rand = RCs(rows, cols);
    case 3
        Rand = RCs(rows, cols);
    otherwise
        Rand = fcn(x0);
end

    function [ Rands ] = RCs(rows, cols)
        Rands = zeros(rows, cols);
        for j = 1 : cols
            for i = 1 : rows
                Rands(i, j) = fcn(x0);
                x0 = feval(OPTIONS_.x0);
            end
        end
    end

end
