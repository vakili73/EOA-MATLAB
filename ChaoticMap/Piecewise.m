function [ Rand ] = Piecewise(rows, cols, OPTIONS)
% where a is a control parameter in (0, 0.5)
a = 0.1;
OPTIONS_.x0 = @rand;
if  nargin == 3 && ~isempty(OPTIONS) && isstruct(OPTIONS)
    if isfield(OPTIONS, 'a'); a = OPTIONS.a; end
    if isfield(OPTIONS, 'x0'); OPTIONS_.x0 = OPTIONS.x0; end
end
x0 = feval(OPTIONS_.x0);

%% Piecewise function
    function [y] = fcn(x0)
        if (0 <= x0) && (x0 <= a)
            y = x0 / a;
            return;
        end
        if (a <= x0) && (x0 <= 0.5)
            y = (x0 - a) / (0.5 - a);
            return;
        end
        if (0.5 <= x0) && (x0 <= 1-a)
            y = (1 - a - x0) / (0.5 - a);
            return;
        end
        if (1-a <= x0) && (x0 <= 1)
            y = (1 - x0) / a;
            return;
        end
        y = fcn(feval(OPTIONS_.x0));
    end

%% random generation method
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
