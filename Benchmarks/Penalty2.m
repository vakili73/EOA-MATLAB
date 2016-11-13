function [InitFunction, CostFunction] = Penalty2
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
x = zeros(1, p);
for popindex = 1 : length(Population)
    Population(popindex).cost = 0;
    for i = 1 : p
        x(i) = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        if (x(i) > 5)
            u = 100 * (x(i) - 5)^4;
        elseif (x(i) < -5)
            u = 100 * (-x(i) - 5)^4;
        else
            u = 0;
        end
        Population(popindex).cost = Population(popindex).cost + u;
    end
    Population(popindex).cost = Population(popindex).cost + 0.1 * ((sin(pi*3*x(1)))^2 + (x(p) - 1)^2 * (1 + (sin(2*pi*x(p)))^2));
    for i = 1 : p-1
        Population(popindex).cost = Population(popindex).cost + 0.1 * ((x(i) - 1)^2 * (1 + (sin(3*pi*x(i+1)))^2));
    end
end
return