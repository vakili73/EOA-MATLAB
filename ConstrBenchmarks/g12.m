function [InitFunction, CostFunction, FeasibleFunction] = g12
% g12 Function 
% Dimension: 3
% Domain: Variable
% min x: 5,5,5
% min f(x): -1
% Number of constraints: 9^3
InitFunction = @g12Init;
CostFunction = @g12Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g12Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -1;
OPTIONS.numVar = 3; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [ 0  0  0];
OPTIONS.MaxDomain = [10 10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g12Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    i = 0;
    Population(popindex).cost = -(100 - ((x(1)-5)^2) - ((x(2)-5)^2) - ((x(3)-5)^2))/100;
    % Compute the level of constraint violation
    % g1 : g729 <= 0
    for p = 1 : 9
        for q = 1 : 9
            for r = 1 : 9
                h = ((x(1) - p)^2) + ((x(2) - q)^2) + ((x(3) - r)^2) - 0.0625;
                if h <= 0
                    Population(popindex).g(1) = h;
                    i = 1;
                end
            end
        end
    end
    if i == 0
        Population(popindex).g(1) = 100;    
    end
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return