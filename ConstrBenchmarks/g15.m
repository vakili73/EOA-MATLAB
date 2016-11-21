function [InitFunction, CostFunction, FeasibleFunction] = g15
% g15 Function 
% Dimensions: 3
% Domain: Variable
% min x: 3.51212..,0.21698..,3.55217.. 
% min f(x): 961.715022...
% Number of Constraints: 2
InitFunction = @g15Init;
CostFunction = @g15Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g15Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 961.715022289961;
OPTIONS.numVar = 3; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [0 0 0];
OPTIONS.MaxDomain = [10 10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
	Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g15Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize 
    x = Population(popindex).chrom;
    Population(popindex).cost = 1000 - x(1)^2 - 2*x(2)^2 - x(3)^2 - x(1)*x(2) - x(1)*x(3);
    % Compute the level of constraint violation
    % h1 = 0
    % Common Method: Convert equality to 2 inequalities
    h1 = x(1)^2 + x(2)^2 + x(3)^2 - 25;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    % h2 = 0
    % Common Method: Convert equality to 2 inequalities
    h2 = 8*x(1) + 14*x(2) + 7*x(3) - 56;
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    % Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return