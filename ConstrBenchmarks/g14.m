function [InitFunction, CostFunction, FeasibleFunction] = g14
% g14 Function
% Dimensions: 10
% Domain: Variable
% min x: 0.04..., 0.14772...,
% min f(x): -47.76... best known
% Constraints: 3
InitFunction = @g14Init;
CostFunction = @g14Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g14Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -47.7;
OPTIONS.numVar = 10; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [ 0  0  0  0  0  0  0  0  0  0];
OPTIONS.MaxDomain = [10 10 10 10 10 10 10 10 10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g14Cost(Population, ~)
popsize = length(Population);
c = [-6.089, -17.164, -34.054, -5.914, -24.721, -14.986, -24.1, -10.708, -26.662, -22.1791];
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    sumx = sum(x);
    Population(popindex).cost = sum(x .* (c + log(x / sumx)));
    % Compute the level of constraint violation
    %Common Method: Convert equality to 2 inequalities
    % h1 = 0
    h1 = x(1) + 2*x(2) + 2*x(3) + x(6) + x(10) - 2;
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    % h2 = 0
    h2 = x(4) + 2*x(5) + x(6) + x(7) - 1;
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;    
    % h3 = 0
    h3 = x(3) + x(7) + x(8) + 2*x(9) + x(10) - 1;
    Population(popindex).g(5) = h3 - epsilon;
    Population(popindex).g(6) = -h3 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return