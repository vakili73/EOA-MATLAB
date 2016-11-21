function [InitFunction, CostFunction, FeasibleFunction] = g03
% g0 Function
% Dimension: variable, n=10
% Domain: 0 <= xi <= 1
% min x: Unknown
% min f(x): 1/sqrt(n)????
% Constraints: 1 nonlinear equality constraint
InitFunction = @g03Init;
CostFunction = @g03Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g03Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue= -1;
OPTIONS.numVar = 10; %Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = zeros(1, OPTIONS.numVar);
OPTIONS.MaxDomain = ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = false;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g03Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
n = OPTIONS.numVar;
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize
    d1 = 1;
    for i = 1 : OPTIONS.numVar
        d1 = d1 * Population(popindex).chrom(i);
    end
    Population(popindex).cost = -(sqrt(n))^n * d1;
    % Compute the level of constraint violation
    % Common Method: Convert equality to 2 inequalities
    h1 = 0;
    for dim = 1 : OPTIONS.numVar
        h1 = h1 + (Population(popindex).chrom(dim))^2;
    end
    h1 = h1 - 1;
    % g1 <= 0
    Population(popindex).g(1) = h1 - epsilon;
    % g2 <= 0
    Population(popindex).g(2) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return