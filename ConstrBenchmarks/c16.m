function [InitFunction, CostFunction, FeasibleFunction] = c16
% c16 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 4
InitFunction = @c16Init;
CostFunction = @c16Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c16Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -10 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +10 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c16Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
epsilon = 0.0001;
for popindex = 1 : popsize
    temp = [0 1 0 1 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    for i = 1 : D
        temp(1) = temp(1) + z(i)^2/4000;
    end
    for i = 1 : D
        temp(2) = temp(2) * cos(z(i)/sqrt(i));
    end
    Population(popindex).cost = temp(1) - temp(2) + 1;
    % Compute the level of constraint violation
    % g1 <= 0
    for i = 1 : D
        temp(3) = temp(3) + z(i)^2 - 100*cos(pi*z(i)) + 10;
    end
    Population(popindex).g(1) = temp(3);
    % g2 <= 0
    for i = 1 : D
        temp(4) = temp(4) * z(i);
    end
    Population(popindex).g(2) = temp(4);
    % h1 = 0
    for i = 1 : D
        temp(5) = temp(5) + z(i)*sin(sqrt(abs(z(i))));
    end
    h1 = temp(5);
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constraints: Bigger the value, bigger the constrain violation
end
return