function [InitFunction, CostFunction] = Griewank
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.MinDomain = -600 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +600 * ones(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true; % the order of the independent variables does matter
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
for popindex = 1 : length(Population)
    sum = 0;
    % Product term fixed on June 1, 2011
    % product = 0;
    product = 1;
    for i = 1 : p
        x = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        sum = sum + x^2;
        product = product * cos(x / sqrt(i));
    end
    Population(popindex).cost = 1 + sum / 4000 - product;
end
return