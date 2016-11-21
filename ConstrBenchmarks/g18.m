function [InitFunction, CostFunction, FeasibleFunction] = g18
% g18 Function 
% Dimension: 9
% Domain: Variable
% min x: -0.657776192427943163,-0.153418773482438542,0.323413871675240938,-0.94625
%          7611651304398,-0.657776194376798906,-0.753213434632691414,0.3234138741235
%          76972,-0.346462947962331735,0.59979466285217542
% min f(x): -0.866025403784439
% Number of constraints: 13
InitFunction = @g18Init;
CostFunction = @g18Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g18Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -0.866025403784439;
OPTIONS.numVar = 9; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [-10 -10 -10 -10 -10 -10 -10 -10  0];
OPTIONS.MaxDomain = [+10 +10 +10 +10 +10 +10 +10 +10 20];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g18Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = -0.5*(x(1)*x(4) - x(2)*x(3) + x(3)*x(9) - x(5)*x(9) + x(5)*x(8) - x(6)*x(7));
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = x(3)^2 + x(4)^2 - 1;
    %g2 <= 0
    Population(popindex).g(2) = x(9)^2 - 1;
    %g3 <= 0
    Population(popindex).g(3) = x(5)^2 + x(6)^2 - 1;
    %g4 <= 0
    Population(popindex).g(4) = x(1)^2 + (x(2) - x(9))^2 - 1;
    %g5 <= 0
    Population(popindex).g(5) = (x(1) - x(5))^2 + (x(2) - x(6))^2 - 1;
    %g6 <= 0
    Population(popindex).g(6) = (x(1) - x(7))^2 + (x(2) - x(8))^2 - 1;
    %g7 <= 0
    Population(popindex).g(7) = (x(3) - x(5))^2 + (x(4) - x(6))^2 - 1;
    %g8 <= 0
    Population(popindex).g(8) = (x(3) - x(7))^2 + (x(4) - x(8))^2 - 1;
    %g9 <= 0
    Population(popindex).g(9) = x(7)^2 + (x(8) - x(9))^2 - 1;
    %g10 <= 0
    Population(popindex).g(10) = x(2)*x(3) - x(1)*x(4);
    %g11 <= 0
    Population(popindex).g(11) = -x(3)*x(9);
    %g12 <= 0
    Population(popindex).g(12) = x(5)*x(9);
    %g13 <= 0
    Population(popindex).g(13) = x(6)*x(7) - x(5)*x(8);
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return