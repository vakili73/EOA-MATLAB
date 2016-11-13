function [InitFunction, CostFunction] = Schwefel222
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
if ~isfield(OPTIONS, 'MinDomain')
    OPTIONS.MinDomain = -10 * ones(1, OPTIONS.numVar);
end
if ~isfield(OPTIONS, 'MaxDomain')
    OPTIONS.MaxDomain = +10 * ones(1, OPTIONS.numVar);
end
OPTIONS.OrderDependent = false; % the order of the independent variables does not matter
OPTIONS = ComputeRandomShift(OPTIONS);
% Initialize population
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
p = length(Population(1).chrom);
if ~exist('OPTIONS', 'var') || isempty(OPTIONS) || ~isfield(OPTIONS, 'ShiftAmount')
    OPTIONS.ShiftAmount = zeros(1, p);
end
for popindex = 1 : length(Population)
    sum = 0;
    product = 1;
    for i = 1 : p
        x = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        sum = sum + abs(x);
        product = product * abs(x);
    end
    Population(popindex).cost = sum + product;
end
return