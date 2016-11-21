function [InitFunction, CostFunction, FeasibleFunction] = g07
% g07 Function 
% Dimensions: 10
% Domain: Variable
% min x: 2.17..., 2.36..., 8.77.., 5.09.., 0.99.., 1.43.., 1.32.., 9.82..,
%       8.28.., 8.37..
% min f(x): 24.306...
% Number of constraints: 8
InitFunction = @g07Init;
CostFunction = @g07Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g07Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 24.306;
OPTIONS.numVar = 10; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = -10 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = 10 * ones(1, OPTIONS.numVar);
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g07Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = x(1)^2 + x(2)^2 + x(1)*x(2) - 14*x(1) - 16*x(2) + (x(3)-10)^2 + 4*(x(4)-5)^2 + ...
        (x(5)-3)^2 + 2*(x(6)-1)^2 + 5*x(7)^2 + 7*(x(8)-11)^2 + 2*(x(9)-10)^2 + (x(10)-7)^2 + 45;
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = -105 + 4*x(1)+5*x(2) - 3*x(7) + 9*x(8);
    %g2 <= 0
    Population(popindex).g(2) = 10*x(1) - 8*x(2) - 17*x(7) + 2*x(8);
    %g3 <= 0
    Population(popindex).g(3) = -8*x(1) + 2*x(2) + 5*x(9) - 2*x(10) - 12;
    %g4 <= 0
    Population(popindex).g(4) = 3*(x(1)-2)^2 + 4*(x(2)-3)^2 + 2*x(3) - 7*x(4) - 120;
    %g5 <= 0
    Population(popindex).g(5) = 5*x(1)^2 + 8*x(2) + (x(3)-6)^2 - 2*x(4) - 40;
    %g6 <= 0
    Population(popindex).g(6) = x(1)^2 + 2*(x(2)-2)^2 - 2*x(1)*x(2) + 14*x(5) - 6*x(6);
    %g7 <= 0
    Population(popindex).g(7) = 0.5*(x(1)-8)^2 + 2*(x(2)-4)^2 + 3*x(5)^2 - x(6) - 30;
    %g8 <= 0
    Population(popindex).g(8) = -3*x(1) + 6*x(2) + 12*(x(9)-8)^2 - 7*x(10);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return