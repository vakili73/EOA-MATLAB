function [InitFunction, CostFunction, FeasibleFunction] = g08
% g08 Function 
% Dimension: 2
% Domain: Variable
% min x: 1.22797...,4.245373......
% min f(x): -0.095825....
% Number of constraints: 2
InitFunction = @g08Init;
CostFunction = @g08Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g08Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -0.095825041;
OPTIONS.numVar = 2; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [0 0];
OPTIONS.MaxDomain = [10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g08Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = -((sin(2*pi*x(1))*sin(2*pi*x(2)))/(x(1)^3*(x(1) + x(2))));
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = x(1)^2 - x(2) + 1;
    %g2 <= 0
    Population(popindex).g(2) = 1 - x(1) + (x(2) - 4)^2;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return