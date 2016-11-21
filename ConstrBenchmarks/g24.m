function [InitFunction, CostFunction, FeasibleFunction] = g24
% g22 Function
% Dimensions: 2
% Domain: Variable
% min x: 2.32952019747762, 3.17849307411774
% min f(x): -5.50801327159536
% Constraints: 2
InitFunction = @g24Init;
CostFunction = @g24Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g24Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -5.50801327159536;
OPTIONS.numVar = 2; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [0 0];
OPTIONS.MaxDomain = [3 4];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g24Cost(Population, ~)
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = -x(1) - x(2);
    % Compute the level of constraint violation
    %g1 <=0
    Population(popindex).g(1) = -2*power(x(1),4) + 8*power(x(1),3) - 8*power(x(1),2) + x(2) - 2;
    %g2 <=0
    Population(popindex).g(2) = -4*power(x(1),4) + 32*power(x(1),3) - 88*power(x(1),2) + 96*x(1) + x(2) - 36;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return