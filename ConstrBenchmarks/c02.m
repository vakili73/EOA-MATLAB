function [InitFunction, CostFunction, FeasibleFunction] = c02
% c02 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 3
InitFunction = @c02Init;
CostFunction = @c02Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c02Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -5.12 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +5.12 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c02Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
epsilon = 0.0001; % inequality epsilon
D = length(Population(1).chrom); % Dimension of the problem
for popindex = 1 : popsize
    temp = [0 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    y = z - 0.5;
    Population(popindex).cost = max(z);
    % Compute the level of constraint violation
    %g1 <= 0
    for i = 1 : D
        temp(1) = temp(1) + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
    Population(popindex).g(1) = 10 - temp(1)/D;
    %g2 <= 0
    Population(popindex).g(2) = temp(1)/D - 15;
    %h1 = 0
    for i = 1 : D
        temp(2) = temp(2) + y(i)^2 - 10*cos(2*pi*y(i)) + 10;
    end
    h1 = temp(2)/D - 20;
    Population(popindex).g(3) = h1 - epsilon;
    Population(popindex).g(4) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return