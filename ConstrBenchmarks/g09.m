function [InitFunction, CostFunction, FeasibleFunction] = g09
% g09 Function 
% Dimension: 7
% Domain: Variable
% min x:2.3304..,1.9513..,-0.4775..,4.3657..,-0.6244..,1.0381..,1.5942..
% min f(x): 680.6300...
% Number of constraints: 4
InitFunction = @g09Init;
CostFunction = @g09Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g09Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 680.63005737;
OPTIONS.numVar = 7; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [-10 -10 -10 -10 -10 -10 -10];
OPTIONS.MaxDomain = [+10 +10 +10 +10 +10 +10 +10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g09Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = (x(1)-10)^2 + 5*((x(2)-12)^2) + x(3)^4 + 3*((x(4)-11)^2)...
        + 10*(x(5)^6) + 7*(x(6)^2) + x(7)^4 - 4*x(6)*x(7) - 10*x(6) - 8*x(7);
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = -127 + 2*(x(1)^2) + 3*(x(2)^4) + x(3) + 4*(x(4)^2) + 5*x(5);
    %g2 <= 0
    Population(popindex).g(2) = -282 + 7*x(1) + 3*x(2) + 10*(x(3)^2) + x(4) - x(5);
    %g3 <= 0
    Population(popindex).g(3) = -196 + 23*x(1) + x(2)^2 + 6*x(6)^2 - 8*x(7);
    %g4 <= 0
    Population(popindex).g(4) = 4*x(1)^2 + x(2)^2 - 3*x(1)*x(2) + 2*x(3)^2 + 5*x(6) - 11*x(7);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return