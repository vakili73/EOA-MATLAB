function [MinCost] = AdaptiveHillClimbing(ProblemFunction, DisplayFlag)
% Adaptive hill climbing algorithm for optimizing a general function.
% INPUTS: (1) ProblemFunction is the handle of the function that returns 
%             the handles of the initialization, cost, and feasibility functions.
%         (2) DisplayFlag says whether or not to display information during iterations and plot results.
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
OPTIONS.popsize = 1; % population size
OPTIONS.MaxFEval = 1000; % function evaluation limit
OPTIONS.pmutate = 0.1; % mutation probability
OPTIONS.numVar = 20; % number of dimensions in cost function
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
feval = OPTIONS.popsize;
fevalArr = inf(size(MinCost));
fevalArr(1) = feval;
GenIndex = 1;
while feval <= OPTIONS.MaxFEval
    Child = Population(1);
    for q = 1 : OPTIONS.numVar
        if rand < OPTIONS.pmutate
            feature = OPTIONS.MinDomain(q) + (OPTIONS.MaxDomain(q) - OPTIONS.MinDomain(q)) * rand;
            Child.chrom(q) = feature;
        end
    end
    Child = OPTIONS.CostFunction(Child, OPTIONS);
    feval = feval + 1;
    if Child.cost < Population(1).cost
        Population(1) = Child;
    end        
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
    GenIndex = GenIndex + 1;
    fevalArr(GenIndex) = feval;
end
indices = find(fevalArr <= OPTIONS.MaxFEval);
MinCost = MinCost(1 : indices(end));
fevalArr = fevalArr(1 : indices(end));
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol, fevalArr);
return