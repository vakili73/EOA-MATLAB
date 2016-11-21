function [InitFunction, CostFunction, FeasibleFunction] = c11
% c11 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 1
InitFunction = @c11Init;
CostFunction = @c11Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c11Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -100 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +100 * ones(1, OPTIONS.numVar);
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
function [Population] = c11Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize
    temp = [0 0 0];
    z = (Population(popindex).chrom - OPTIONS.ShiftAmount) * OPTIONS.RotationMatrix;
    y = Population(popindex).chrom - OPTIONS.ShiftAmount + 1;
    for i = 1 : D
        temp(1) = temp(1) - z(i)*cos(2*sqrt(abs(z(i))));
    end
    Population(popindex).cost = temp(1) / D;
    % Compute the level of constraint violation
    % h1 = 0
    for i = 1 : (D-1)
        temp(2) = temp(2) + 100*((y(i)^2 - y(i+1))^2) + (y(i) - 1)^2;
    end
    h1 = temp(2);
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return
