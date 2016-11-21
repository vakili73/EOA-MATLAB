function [InitFunction, CostFunction, FeasibleFunction] = g02
% g02 Function
% Dimension: variable, n=20
% Domain: 0 <= xi <= 10
% min x: Unknown
% min f(x): Unknown, best found: -0.803619
% Constrains: 2: 1 lin ineq and 1 nl ineq
InitFunction = @g02Init;
CostFunction = @g02Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g02Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue= -0.803619;
OPTIONS.numVar = 20; %Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = zeros(1, OPTIONS.numVar);
OPTIONS.MaxDomain = 10 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g02Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    d1 = 0;
    d2 = 1;
    d3 = 0;
    for i = 1 : OPTIONS.numVar
        gene = Population(popindex).chrom(i);
        d1 = d1 + (cos(gene))^4;
        d2 = d2 * (cos(gene))^2;
        d3 = d3 + i * gene * gene;
    end
    Population(popindex).cost = -abs((d1 - 2 * d2) / sqrt(d3));
    % Compute the level of constraint violation
    d1 = 1;
    d2 = 0;
    for dim = 1 : OPTIONS.numVar
        gene = Population(popindex).chrom(dim);
        d1 = d1 * gene;
        d2 = d2 + gene;
    end
    %g1 <= 0
    Population(popindex).g(1) = 0.75 - d1;
    %g2 <= 0
    Population(popindex).g(2) = d2 - 7.5 * OPTIONS.numVar;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return