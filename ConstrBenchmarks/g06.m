function [InitFunction, CostFunction, FeasibleFunction] = g06
% g06 Function 
% Dimensions: 2
% Domain: Variable
% min x: 14.095,0.84296
% min f(x): -6961.81388
% Number of Constraints: 2
InitFunction = @g06Init;
CostFunction = @g06Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g06Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -6961.81388;
OPTIONS.numVar = 2; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [13 0];
OPTIONS.MaxDomain = [100 100];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
	Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g06Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize   
    Population(popindex).cost = (Population(popindex).chrom(1) - 10)^3 + ...
        (Population(popindex).chrom(2) - 20)^3;
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = -(Population(popindex).chrom(1) - 5)^2 - ...
        (Population(popindex).chrom(2) - 5)^2 + 100;
    %g2 <= 0
    Population(popindex).g(2) = (Population(popindex).chrom(1) - 6)^2 + ...
        (Population(popindex).chrom(2) - 5)^2 - 82.81;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return