function [InitFunction, CostFunction] = Penalty1
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.MinDomain = -50 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +50 * ones(1, OPTIONS.numVar);
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
y = zeros(1, p);
for popindex = 1 : length(Population)
    Population(popindex).cost = 0;
    for i = 1 : p
        x = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        y(i) = 1 + (x + 1) / 4;
        if (x > 10)
            u = 100 * (x - 10)^4;
        elseif (x < -10)
            u = 100 * (-x - 10)^4;
        else
            u = 0;
        end
        Population(popindex).cost = Population(popindex).cost + u;
    end
    Population(popindex).cost = Population(popindex).cost + (10 * (sin(pi*y(1)))^2 + (y(p) - 1)^2) * pi / 30;
    for i = 1 : p-1
        Population(popindex).cost = Population(popindex).cost + ((y(i) - 1)^2 * (1 + 10 * (sin(pi*y(i+1)))^2)) * pi / 30;
    end
end
return