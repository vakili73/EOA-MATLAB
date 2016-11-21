function [InitFunction, CostFunction, FeasibleFunction] = c08
% c08 Function 
% Dimension: user defined
% Domain: Variable
% min x: ?
% min f(x): ?
% Number of constraints: 1
InitFunction = @c08Init;
CostFunction = @c08Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = c08Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -15;
% Boundaries for each variable
OPTIONS.MinDomain = -140 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +140 * ones(1, OPTIONS.numVar);
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
function [Population] = c08Cost(Population, OPTIONS)
% Compute the cost of each member in Population
popsize = length(Population);
D = length(Population(1).chrom); % Dimension of the problem
for popindex = 1 : popsize
    temp = [0 0 0];
    z = Population(popindex).chrom - OPTIONS.ShiftAmount + 1;
    y = (Population(popindex).chrom - OPTIONS.ShiftAmount) * OPTIONS.RotationMatrix;
    for i = 1 : (D-1)
        temp(1) = temp(1) + 100*((z(i)^2 - z(i+1))^2) + (z(i) - 1)^2;
    end
    Population(popindex).cost = temp(1);
    % Compute the level of constraint violation
    % g1 = 0
    for i = 1 : D
        temp(2) = temp(2) + y(i)^2;
    end
    for i = 1 : D
        temp(3) = temp(3) + cos(0.1*y(i));
    end
    Population(popindex).g(1) = 0.5 - exp(-0.1*(sqrt(temp(2)/D))) - 3*exp(temp(3)/D) + exp(1);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constraint violation
end
return