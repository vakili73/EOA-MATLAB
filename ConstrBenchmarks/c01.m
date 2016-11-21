function [InitFunction, CostFunction, FeasibleFunction] = c01
% c01 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 2
InitFunction = @c01Init;
CostFunction = @c01Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c01Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = zeros(1, OPTIONS.numVar);
OPTIONS.MaxDomain = 10 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c01Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
for popindex = 1 : popsize
    temp = [0 1 0 1 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    for i = 1 : D
        temp(1) = temp(1) + cos(z(i))^4;
    end
    for i = 1 : D
        temp(2) = temp(2) * cos(z(i))^2;
    end
    for i = 1 : D
        temp(3) = temp(3) + i * z(i)^2;
    end
    Population(popindex).cost = -abs((temp(1) - 2*temp(2))/sqrt(temp(3)));
    % Compute the level of constraint violation
    %g1 <= 0
    for i = 1 : D
        temp(4) = temp(4) * z(i);
    end
    Population(popindex).g(1) = 0.75 - temp(4);
    %g2 <= 0
    for i = 1 : D
        temp(5) = temp(5) + z(i);
    end
    Population(popindex).g(2) = temp(5) - 7.5*D;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return