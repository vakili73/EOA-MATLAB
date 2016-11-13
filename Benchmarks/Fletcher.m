function [InitFunction, CostFunction] = Fletcher
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
global a b alph A
OPTIONS.MinDomain = -pi * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +pi * ones(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true; % the order of the independent variables matters
OPTIONS = ComputeRandomShift(OPTIONS);
a = 200 * rand(OPTIONS.numVar,OPTIONS.numVar) - 100;
b = 200 * rand(OPTIONS.numVar,OPTIONS.numVar) - 100;
alph = 2 * pi * rand(OPTIONS.numVar,1) - pi;
A = zeros(OPTIONS.numVar, 1);
for i = 1 : OPTIONS.numVar
    for j = 1 : OPTIONS.numVar
        A(i) = A(i) + a(i,j) * sin(alph(j)) + b(i,j) * cos(alph(j));
    end
end
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
global a b A
p = length(Population(1).chrom);
for popindex = 1 : length(Population)
    x = zeros(p, 1);
    for j = 1 : p
        x(j) = Population(popindex).chrom(j) - OPTIONS.ShiftAmount(j);
    end
    B = zeros(p, 1);
    for i = 1 : p
        for j = 1 : p
            B(i) = B(i) + a(i,j) * sin(x(j)) + b(i,j) * cos(x(j));
        end
    end
    Population(popindex).cost = 0;
    for i = 1 : p
        Population(popindex).cost = Population(popindex).cost + (A(i) - B(i))^2;
    end
end
return