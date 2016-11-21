function [InitFunction, CostFunction, FeasibleFunction] = c06
% c06 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 2
InitFunction = @c06Init;
CostFunction = @c06Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c06Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -500 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +500 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
OPTIONS.OrderDependent = true;
OPTIONS = ComputeRandomShift(OPTIONS);
if OPTIONS.RotationFlag > 0
    OPTIONS.RotationMatrix = createRotMatrix(OPTIONS.numVar);
else
    OPTIONS.RotationMatrix = eye(OPTIONS.numVar);
end
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = c06Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize
    temp = [0 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount;
    y = (z + 483.6106156535) * OPTIONS.RotationMatrix - 483.6106156535;
    Population(popindex).cost = max(z);
    % Compute the level of constraint violation
    % h1 = 0
    for i = 1 : D
        temp(1) = temp(1) - y(i)*sin(sqrt(abs(y(i))));
    end
    h1 = temp(1) / D;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    % h2 = 0
    for i = 1 : D
        temp(2) = temp(2) - y(i)*cos(0.5*sqrt(abs(y(i))));
    end
    h2 = temp(2) / D;
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;    
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return