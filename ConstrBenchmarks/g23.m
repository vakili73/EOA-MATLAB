function [InitFunction, CostFunction, FeasibleFunction] = g23
% g22 Function
% Dimensions: 9
% Domain: Variable
% min x: 0.00510000000000259465, 99.9947000000000514,
% 9.01920162996045897e-18, 99.9999000000000535, 0.000100000000027086086,
% 2.75700683389584542e-14, 99.9999999999999574, 2000.0100000100000100008
% min f(x): -400.055099999999584
% Constraints: 6
InitFunction = @g23Init;
CostFunction = @g23Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g23Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -400.055099999999584;
OPTIONS.numVar = 9; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [  0   0   0   0   0   0   0   0 0.001];
OPTIONS.MaxDomain = [300 300 100 200 100 300 100 200 0.03];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g23Cost(Population, ~)
popsize = length(Population);
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = -9*x(5) - 15*x(8) + 6*x(1) + 16*x(2) + 10*(x(6) + x(7));
    % Compute the level of constraint violation
    %g1 <=0
    Population(popindex).g(1) = x(9)*x(3) + 0.02*x(6) - 0.025*x(5);
    %g2 <=0
    Population(popindex).g(2) = x(9)*x(4) + 0.02*x(7) - 0.015*x(8);
    %Common Method: Convert equality to 2 inequalities
    % h1 = 0
    h1 = x(1) + x(2) - x(3) - x(4);
    Population(popindex).g(3) = h1 - epsilon;
    Population(popindex).g(4) = -h1 - epsilon;
    % h2 = 0
    h2 = 0.03*x(1) + 0.01*x(2) - x(9)*(x(3) + x(4));
    Population(popindex).g(5) = h2 - epsilon;
    Population(popindex).g(6) = -h2 - epsilon;    
    % h3 = 0
    h3 = x(3) + x(6) - x(5);
    Population(popindex).g(7) = h3 - epsilon;
    Population(popindex).g(8) = -h3 - epsilon;
    % h4 = 0
    h4 = x(4) + x(7) - x(8);
    Population(popindex).g(9) = h4 - epsilon;
    Population(popindex).g(10) = -h4 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return