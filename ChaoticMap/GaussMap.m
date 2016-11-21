function [ Rand ] = GaussMap(rows, cols, OPTIONS)
% Gauss map function is defined by following equation
OPTIONS_.x0 = @rand;
a = 8.75;
b = -0.58;
if  nargin == 3 && ~isempty(OPTIONS) && isstruct(OPTIONS)
    if isfield(OPTIONS, 'x0'); OPTIONS_.x0 = OPTIONS.x0; end
    if isfield(OPTIONS, 'a'); a = OPTIONS.a; end
    if isfield(OPTIONS, 'b'); b = OPTIONS.b; end
end
x0 = feval(OPTIONS_.x0);

%% Gauss function
fcn = @ (x0) exp(-a * (x0^2)) + b;

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
