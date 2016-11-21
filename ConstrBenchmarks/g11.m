function [InitFunction, CostFunction, FeasibleFunction] = g11
% g11 Function 
% Dimensiona: 2
% Domain: Variable
% min x: +- 1/sqrt(2), 1/2
% min f(x): 0.75
% Number of Constraints: 1
InitFunction = @g11Init;
CostFunction = @g11Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g11Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 0.75;
OPTIONS.numVar = 2; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [-1 -1];
OPTIONS.MaxDomain = [1 1];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
	Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g11Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize   
    Population(popindex).cost = (Population(popindex).chrom(1))^2 + (Population(popindex).chrom(2) - 1)^2;
    % Compute the level of constraint violation
    % h1 = 0
    % Common Method: Convert equality to 2 inequalities
    h1 = Population(popindex).chrom(2) - (Population(popindex).chrom(1))^2;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    % Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return