function [InitFunction, CostFunction, FeasibleFunction] = AckleyScaled
InitFunction = @Init;
CostFunction = @Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
% Initialize the population
OPTIONS.MinDomain = -30;
OPTIONS.MaxDomain = +30;
OPTIONS.OrderDependent = false; % the order of the independent variables does not matter
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) * rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, ~)
% Compute the cost of each member in Population
p = length(Population(1).chrom);
for popindex = 1 : length(Population)
    sum1 = 0;
    sum2 = 0;
    for i = 1 : p
        x = Population(popindex).chrom(i);
        % Make the function less sensitive to the even-numbered parameters
        if mod(i, 2) == 0
            x = x / 4;
        end
        sum1 = sum1 + x^2;
        sum2 = sum2 + cos(2*pi*x);
    end
    Population(popindex).cost = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/p)) - exp(sum2/p);    
end
return