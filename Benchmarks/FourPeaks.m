function [InitFunction, CostFunction] = FourPeaks
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.Dim = 40;
OPTIONS.numVar = OPTIONS.Dim;
OPTIONS.OrderDependent = true;
OPTIONS.T = round(OPTIONS.Dim / 10);
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
    head = 0; % number of consecutive 1's at the start of the bit string
    for i = 1 : OPTIONS.Dim
        if Population(popindex).chrom(i) == 1
            head = head + 1;
        else
            break
        end
    end
    tail = 0; % number of consecutive 0's at the end of the bit string
    for i = OPTIONS.Dim : -1 : 1
        if Population(popindex).chrom(i) == 0
            tail = tail + 1;
        else
            break
        end
    end
    fitness = OPTIONS.Dim * ((tail > OPTIONS.T) && (head > OPTIONS.T)) + max(head, tail);
    Population(popindex).cost = 2 * OPTIONS.Dim - fitness;
end
return