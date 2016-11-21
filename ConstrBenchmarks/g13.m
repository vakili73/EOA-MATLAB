function [InitFunction, CostFunction, FeasibleFunction] = g13
% g11 Function 
% Dimensions: 5
% Domain: Variable
% min x: - 1.7171..,1.5957..,1.8272..,-0.7636..,-0.7636.. 
% min f(x): 0.05394...
% Number of Constraints: 3
InitFunction = @g13Init;
CostFunction = @g13Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g13Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 0.053941514;
OPTIONS.numVar = 5; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [-2.3 -2.3 -3.2 -3.2 -3.2];
OPTIONS.MaxDomain = [+2.3 +2.3 +3.2 +3.2 +3.2];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
	Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g13Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
epsilon = 0.0001; % inequality epsilon
for popindex = 1 : popsize 
    x = Population(popindex).chrom;
    Population(popindex).cost = exp(x(1)*x(2)*x(3)*x(4)*x(5));
    % Compute the level of constraint violation
    % h1 = 0
    % Common Method: Convert equality to 2 inequalities
    h1 = x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 - 10;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    % h2 = 0
    % Common Method: Convert equality to 2 inequalities
    h2 = x(2)*x(3) - 5*x(4)*x(5);
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;
    % h3 = 0
    % Common Method: Convert equality to 2 inequalities
    h3 = x(1)^3 + x(2)^3 + 1;
    Population(popindex).g(5) = h3 - epsilon;
    Population(popindex).g(6) = -h3 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    % Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return