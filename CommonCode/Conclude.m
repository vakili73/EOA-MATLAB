function Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol, xAxisArray)
% Output results of evolutionary optimization algorithm.
if DisplayFlag
    % Count the number of unique individuals
    NumUnique = 0;
    DupIndices = [];
    for i = 1 : OPTIONS.popsize
        if find(DupIndices == i), continue, end
        Chrom1 = Population(i).chrom;
        if ~OPTIONS.OrderDependent
            Chrom1 = sort(Chrom1);
        end
        UniqueFlag = 1;
        for j = i+1 : OPTIONS.popsize
            Chrom2 = Population(j).chrom;
            if ~OPTIONS.OrderDependent
                Chrom2 = sort(Chrom2);
            end
            if isequal(Chrom1, Chrom2)
                DupIndices = [DupIndices j]; %#ok<AGROW>
                UniqueFlag = 0;
                break
            end
        end
        NumUnique = NumUnique + UniqueFlag;
    end
    disp([num2str(NumUnique), ' unique individuals in final population.']);
    % Display the best solution
    Chrom = Population(1).chrom;
    if ~OPTIONS.OrderDependent
        Chrom = sort(Chrom);
    end
    disp(['Best chromosome = ', num2str(Chrom)]); 
    % Plot the minimum and average cost.
    % For constrained cost functions, only plot the cost values for the generation numbers after 
    % which the EA found feasible solutions.
    if ~exist('xAxisArray', 'var') || isempty(xAxisArray)
        xAxisArray = 0 : length(MinCost)-1;
        xlabelString = 'Generation';
    else
        xlabelString = 'Function Evaluations';
    end
    FirstFeasibleIndex = find(MinConstrViol == 0, 1, 'first');
    if ~isempty(FirstFeasibleIndex)
        figure, hold on
        if OPTIONS.popsize > 1
            plot(xAxisArray(FirstFeasibleIndex:end), AvgCost(FirstFeasibleIndex:end), 'r--', ...
                xAxisArray(FirstFeasibleIndex:end), MinCost(FirstFeasibleIndex:end), 'b-')
            legend('Average Cost', 'Minimum Cost')
        else
            plot(xAxisArray(FirstFeasibleIndex:end), MinCost(FirstFeasibleIndex:end), 'b-')
        end
        tempAxis = axis;
        axis([0 tempAxis(2:4)])
        xlabel(xlabelString)
        ylabel('Cost')
    end
    % If the cost function is constrained, plot the number of constraint violations
    if isfield(Population, 'G')
        figure, hold on
        if OPTIONS.popsize > 1
            plot(xAxisArray, AvgConstrViol, 'r--', xAxisArray, MinConstrViol, 'b-')
            legend('Average # Constr. Violations', 'Minimum # Constr. Violations')
        else
            plot(xAxisArray, MinConstrViol, 'b-')
        end
        xlabel(xlabelString)
        ylabel('Number of Constraint Violations')
    end
end
return