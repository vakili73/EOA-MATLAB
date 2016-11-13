function [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ...
    ComputeCostAndConstrViol(Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag)
% Compute the minimum cost, the average cost, the minimum number of constraint violations, and
% the average number of constraint violations in the population
if isfield(Population, 'G')
    % If the problem is constrained, the minimum and average costs are taken over feasible solutions
    Feasibles = Population([Population.G] == 0);
    if isempty(Feasibles)
        MinCost(GenIndex+1) = NaN;
    else
        MinCost(GenIndex+1) = Feasibles(1).cost;
    end
    AvgCost(GenIndex+1) = mean([Feasibles.cost]);
    popSize = length(Population);
    numConstr = length(Population(1).g);
    ConstrViolAmt = reshape([Population.g], numConstr, popSize)';
    ConstrViolFlags = (ConstrViolAmt ~= 0);
    NumConstrViol = sum(ConstrViolFlags, 2);
    MinConstrViol(GenIndex+1) = min(NumConstrViol);
    AvgConstrViol(GenIndex+1) = mean(NumConstrViol);
else
    MinCost(GenIndex+1) = Population(1).cost;
    AvgCost(GenIndex+1) = mean([Population.cost]);
end
if DisplayFlag
    % Display info to screen
    if isfield(Population, 'G')
        disp(['Gen # ', num2str(GenIndex), ': min cost = ', ...
            num2str(MinCost(GenIndex+1)), ', ave cost = ', num2str(AvgCost(GenIndex+1)), ...
            ', min # viol = ', num2str(MinConstrViol(GenIndex+1)), ...
            ', ave # viol = ', num2str(AvgConstrViol(GenIndex+1))]);
    else
        disp(['Gen # ', num2str(GenIndex), ': min cost = ', ...
            num2str(MinCost(GenIndex+1)), ', ave cost = ', num2str(AvgCost(GenIndex+1))]);
    end
end