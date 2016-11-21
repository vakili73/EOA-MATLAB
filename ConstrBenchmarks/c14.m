function [InitFunction, CostFunction, FeasibleFunction] = c14
% c14 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 3
InitFunction = @c14Init;
CostFunction = @c14Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c14Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -1000 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +1000 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c14Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
for popindex = 1 : popsize
    temp = [0 0 0];
    y = Population(popindex).chrom - OPTIONS.ShiftAmount;
    z = y + 1;
    for i = 1 : (D-1)
        temp(1) = temp(1) + (100*(z(i)^2-z(i+1))^2) + (z(i)-1)^2;
    end
    Population(popindex).cost = temp(1);
    % Compute the level of constraint violation
    % g1 <= 0
    for i = 1 : D
        temp(2) = temp(2) + y(i)*cos(sqrt(abs(y(i))));
    end
    Population(popindex).g(1) = -temp(2) - D;
    % g2 <= 0
    Population(popindex).g(2) = temp(2) - D;
    % g3 <= 0
    for i = 1 : D
        temp(3) = temp(3) + y(i)*sin(sqrt(abs(y(i))));
    end
    Population(popindex).g(3) = temp(3) - 10*D;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constraints: Bigger the value, bigger the constrain violation
end
return
