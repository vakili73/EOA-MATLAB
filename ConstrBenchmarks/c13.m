function [InitFunction, CostFunction, FeasibleFunction] = c13
% c13 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 3
InitFunction = @c13Init;
CostFunction = @c13Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c13Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -500 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +500 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c13Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
for popindex = 1 : popsize
    temp = [0 0 0 0 1];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    for i = 1 : D
        temp(1) = temp(1) - z(i)*sin(sqrt(abs(z(i))));
    end
    Population(popindex).cost = temp(1);
    % Compute the level of constraint violation
    % g1 <= 0
    for i = 1 : D
        temp(2) = temp(2) + z(i)^2;
    end
    Population(popindex).g(1) = -50 + temp(2) / 100 / D;
    % g2 <= 0
    for i = 1 : D
        temp(3) = temp(3) + sin(pi*z(i)/50);
    end
    Population(popindex).g(2) = 50 * temp(3) / D;
    % g3 <= 0
    for i = 1 : D
        temp(4) = temp(4) + z(i)^2/4000;
    end
    for i = 1 : D
        temp(5) = temp(5) * cos(z(i)/sqrt(i));
    end
    Population(popindex).g(3) = 75 - 50*(temp(4) - temp(5) + 1);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return