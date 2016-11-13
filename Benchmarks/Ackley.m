function [InitFunction, CostFunction] = Ackley
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
% Initialize the population
if ~isfield(OPTIONS, 'MinDomain')
    OPTIONS.MinDomain = -30 * ones(1, OPTIONS.numVar);
end
if ~isfield(OPTIONS, 'MaxDomain')
    OPTIONS.MaxDomain = +30 * ones(1, OPTIONS.numVar);
end
OPTIONS.OrderDependent = false; % the order of the independent variables does not matter
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
p = length(Population(1).chrom);
for popindex = 1 : length(Population)
    sum1 = 0;
    sum2 = 0;
    for i = 1 : p
        x = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        sum1 = sum1 + x^2;
        sum2 = sum2 + cos(2*pi*x);
    end
    Population(popindex).cost = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/p)) - exp(sum2/p);    
end
return