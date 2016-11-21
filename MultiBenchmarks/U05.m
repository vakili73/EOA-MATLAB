function [InitFunction, CostFunction, FeasibleFunction] = U05
% Adapated from cec09.m, which was downloaded from www.ntu.edu.sg/home/EPNSugan/
InitFunction = @Init;
CostFunction = @Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
% Initialize the population
if ~isfield(OPTIONS, 'MinDomain')
    OPTIONS.MinDomain = [0, -1 * ones(1, OPTIONS.numVar-1)];
end
if ~isfield(OPTIONS, 'MaxDomain')
    OPTIONS.MaxDomain = ones(1, OPTIONS.numVar);
end
OPTIONS.OrderDependent = true; % the order of the independent variables does matter
OPTIONS = ComputeRandomShift(OPTIONS);
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
num = length(Population);
dim =  OPTIONS.numVar;
y = zeros(1, dim);
for pop = 1 : num
    x = Population(pop).chrom;
    epsi = 0.1;
    N = 10;
    for j = 2 : 1 : dim
        y(j) = x(j) - sin((6*pi*x(1)) + (j*pi/dim));
    end
    % f1
    temp1 = 0;
    for j = 3 : 2 : dim
        temp1 = temp1 + (2*(y(j)^2)) - cos(4*pi*y(j)) + 1;
    end
    J1 = 3 : 2 : dim;
    Population(pop).cost(1) = x(1) + (((1/2*N) + epsi)*abs(sin(2*N*pi*x(1)))) + (2 * temp1 / length(J1));
    % f2
    temp2 = 0;
    for j = 2 : 2 : dim
        temp2 = temp2 + (2*(y(j)^2)) - cos(4*pi*y(j)) + 1;
    end
    J2 = 2 : 2 : dim;
    Population(pop).cost(2) = 1 - x(1) + (((1/2*N) + epsi)*abs(sin(2*N*pi*x(1)))) + (2 * temp2 / length(J2));
end
return