function [InitFunction, CostFunction] = Schwefel12
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.MinDomain = -65.536 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +65.536 * ones(1, OPTIONS.numVar);
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
    Population(popindex).cost = 0;
    for i = 1 : p
        innersum = 0;
        for j = 1 : i
            innersum = innersum + Population(popindex).chrom(j) - OPTIONS.ShiftAmount(j);
        end
        Population(popindex).cost = Population(popindex).cost + innersum^2;
    end
end
return