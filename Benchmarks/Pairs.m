function [InitFunction, CostFunction] = Pairs
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.Dim = 40; % must be even
OPTIONS.numVar = OPTIONS.Dim;
OPTIONS.OrderDependent = true;
% Initialize population
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    Population(popindex).chrom = round(rand(1, OPTIONS.numVar));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(OPTIONS, Population)
% Compute the cost of each member in Population
for popindex = 1 : length(Population)
    pairs = 0; % number of alternating pairs of 0's and pairs of 1's
    incremented = true;
    for i = 1 : OPTIONS.Dim / 2
        if Population(popindex).chrom(2*i-1) == Population(popindex).chrom(2*i)
            if (i == 1) || ~incremented
                pairs = pairs + 1;
                incremented = true;
            elseif Population(popindex).chrom(2*i) ~= Population(popindex).chrom(2*i-2)
                pairs = pairs + 1;
                incremented = true;
            else
                incremented = false;
            end
        else
            incremented = false;
        end
    end
    Population(popindex).cost = OPTIONS.Dim / 2 - pairs;
end
return