function [InitFunction, CostFunction, FeasibleFunction] = c10
% c10 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 1
InitFunction = @c10Init;
CostFunction = @c10Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c10Init(OPTIONS)
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
function [Population] = c10Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
epsilon = 0.0001;
for popindex = 1 : popsize
    temp = [0 0];
    y = Population(popindex).chrom - OPTIONS.ShiftAmount + 1;
    z = (Population(popindex).chrom - OPTIONS.ShiftAmount) * OPTIONS.RotationMatrix;
    for i = 1 : (D-1)
        temp(1) = temp(1) + 100*(z(i)^2 - z(i+1))^2 + (z(i)-1)^2;
    end
    Population(popindex).cost = temp(1);
    % Compute the level of constraint violation
    % h1 = 0
    for i = 1 : D
        temp(2) = temp(2) + y(i)*sin(sqrt(abs(y(i))));
    end
    h1 = temp(2);
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    % Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    % Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return
 