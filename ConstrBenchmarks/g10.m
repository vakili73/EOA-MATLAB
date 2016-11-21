function [InitFunction, CostFunction, FeasibleFunction] = g10
% g10 Function 
% Dimension: 8
% Domain: Variable
% min
% x:579.3066..,1359.9706..,5109.9706..,182.0176..,295.6011..,217.9823..,
%   286.4165..,385.6011..
% min f(x): 7049.2480...
% Number of constraints: 6
InitFunction = @g10Init;
CostFunction = @g10Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g10Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 7049.24802052;
OPTIONS.numVar = 8; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [  100  1000  1000   10   10   10   10   10];
OPTIONS.MaxDomain = [10000 10000 10000 1000 1000 1000 1000 1000];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g10Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = x(1) + x(2) + x(3);
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = -1 + 0.0025*(x(4) + x(6));
    %g2 <= 0
    Population(popindex).g(2) = -1 + 0.0025*(x(5)+x(7)-x(4));
    %g3 <= 0
    Population(popindex).g(3) = -1 + 0.01*(x(8) - x(5));
    %g4 <= 0
    Population(popindex).g(4) = -x(1)*x(6) + 833.33252*x(4) + 100*x(1) - 83333.333;
    %g5 <= 0
    Population(popindex).g(5) = -x(2)*x(7) + 1250*x(5) + x(2)*x(4) - 1250*x(4);
    %g6 <= 0
    Population(popindex).g(6) = -x(3)*x(8) + 1250000 + x(3)*x(5) - 2500*x(5);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return