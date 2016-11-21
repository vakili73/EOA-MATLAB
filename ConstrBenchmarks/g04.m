function [InitFunction, CostFunction, FeasibleFunction] = g04
% g04 Function 
% Dimension: 5
% Domain: Variable
% min x: 78,33,29.995 256 025 682, 45, 36.775 812 905 788
% min f(x): -30665.539
% Number of  Constraints: 6
InitFunction = @g04Init;
CostFunction = @g04Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g04Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue= -30665.539;
OPTIONS.numVar = 5; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [78 33 27 27 27];
OPTIONS.MaxDomain = [102 45 45 45 45];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
	Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g04Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize   
    Population(popindex).cost = 5.3578547 * (Population(popindex).chrom(3))^2 + ...
        0.8356891 * Population(popindex).chrom(1) * Population(popindex).chrom(5) + ...
        37.293239 * Population(popindex).chrom(1) - 40792.141;
    % Compute the level of constraint violation
    % g1 <= 0
    Population(popindex).g(1) = 85.334407 + 0.0056858*Population(popindex).chrom(2)*Population(popindex).chrom(5) + ...
        0.0006262*Population(popindex).chrom(1)*Population(popindex).chrom(4) -...
        0.0022053*Population(popindex).chrom(3)*Population(popindex).chrom(5) - 92;
    % g2 <= 0
    Population(popindex).g(2) = -85.334407 - 0.0056858*Population(popindex).chrom(2)*Population(popindex).chrom(5) - ...
        0.0006262*Population(popindex).chrom(1)*Population(popindex).chrom(4) -...
        0.0022053*Population(popindex).chrom(3)*Population(popindex).chrom(5);
    % g3 <= 0
    Population(popindex).g(3) = 80.51249 + 0.0071317*Population(popindex).chrom(2)*Population(popindex).chrom(5) + ...
        0.0029955*Population(popindex).chrom(1)*Population(popindex).chrom(2) + ...
        0.0021813*Population(popindex).chrom(3)*Population(popindex).chrom(3) - 110;
    % g4 <= 0
    Population(popindex).g(4) = -80.51249 - 0.0071317*Population(popindex).chrom(2)*Population(popindex).chrom(5) - ...
        0.0029955*Population(popindex).chrom(1)*Population(popindex).chrom(2) - ...
        0.0021813*Population(popindex).chrom(3)*Population(popindex).chrom(3) + 90;
    % g5 <= 0 
    Population(popindex).g(5) = 9.300961 + 0.0047026*Population(popindex).chrom(3)*Population(popindex).chrom(5) + ...
        0.0012547*Population(popindex).chrom(1)*Population(popindex).chrom(3) + ...
        0.0019085*Population(popindex).chrom(3)*Population(popindex).chrom(4) - 25;
    % g6 <= 0
    Population(popindex).g(6) = -9.300961 - 0.0047026*Population(popindex).chrom(3)*Population(popindex).chrom(5) - ...
        0.0012547*Population(popindex).chrom(1)*Population(popindex).chrom(3) - ...
        0.0019085*Population(popindex).chrom(3)*Population(popindex).chrom(4) + 20;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return