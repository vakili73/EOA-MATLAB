function [InitFunction, CostFunction, FeasibleFunction] = c05
% c05 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 2
InitFunction = @c05Init;
CostFunction = @c05Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c05Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -600 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +600 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c05Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    temp = [0 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    Population(popindex).cost = max(z);
    % Compute the level of constraint violation
    %h1 = 0
    for i = 1 : D
        temp(1) = temp(1) - z(i)*sin(sqrt(abs(z(i))));
    end
    h1 = temp(1) / D;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    %h2 = 0
    for i = 1 : D
        temp(2) = temp(2) - z(i)*cos(0.5*sqrt(abs(z(i))));
    end
    h2 = (1/D)*temp(2);
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;    
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return