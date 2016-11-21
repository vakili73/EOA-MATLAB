function [InitFunction, CostFunction, FeasibleFunction] = g01
% g01 Function 
% Dimension: 13
% Domain: Variable
% min x: 1,1,1,1,1,1,1,1,1,3,3,3,1
% min f(x): -15
% Number of constraints: 9
InitFunction = @g01Init;
CostFunction = @g01Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g01Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
OPTIONS.numVar = 13; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [0 0 0 0 0 0 0 0 0 0 0 0 0];
OPTIONS.MaxDomain = [1 1 1 1 1 1 1 1 1 100 100 100 1];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g01Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    d1 = 0;
    d2 = 0;
    for i = 1 : 4
        gene = Population(popindex).chrom(i);
        d1 = d1 + gene;
        d2 = d2 + gene^2;
    end
    p1 = 5 * d1;
    p2 = 5 * d2;
    d3 = 0;
    for i = 5 : 13
        gene = Population(popindex).chrom(i);
        d3 = d3 + gene;
    end
    p3 = d3;
    Population(popindex).cost = p1 - p2 - p3;
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = 2*Population(popindex).chrom(1) + 2*Population(popindex).chrom(2) + ...
        Population(popindex).chrom(10) + Population(popindex).chrom(11) - 10;
    %g2 <= 0
    Population(popindex).g(2) = 2*Population(popindex).chrom(1) + 2*Population(popindex).chrom(3) + ...
        Population(popindex).chrom(10) + Population(popindex).chrom(12) - 10;
    %g3 <= 0
    Population(popindex).g(3) = 2*Population(popindex).chrom(2) + 2*Population(popindex).chrom(3) + ...
        Population(popindex).chrom(11) + Population(popindex).chrom(12) - 10;
    %g4 <= 0
    Population(popindex).g(4) = -8*Population(popindex).chrom(1) + Population(popindex).chrom(10);
    %g5 <= 0
    Population(popindex).g(5) = -8*Population(popindex).chrom(2) + Population(popindex).chrom(11);
    %g6 <= 0
    Population(popindex).g(6) = -8*Population(popindex).chrom(3) + Population(popindex).chrom(12);
    %g7 <= 0
    Population(popindex).g(7) = -2*Population(popindex).chrom(4) - Population(popindex).chrom(5) + ...
        Population(popindex).chrom(10);
    %g8 <= 0
    Population(popindex).g(8) = -2*Population(popindex).chrom(6) - Population(popindex).chrom(7) + ...
        Population(popindex).chrom(11);
    %g9 <= 0
    Population(popindex).g(9) = -2*Population(popindex).chrom(8) - Population(popindex).chrom(9) + ...
        Population(popindex).chrom(12);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return