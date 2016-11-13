function [MinCost] = NextHillClimbing(ProblemFunction, DisplayFlag)
% Next ascent hill climbing algorithm for optimizing a general function.
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
Children(OPTIONS.numVar) = struct(Population);
Best = Population(1);
GenIndex = 1;
while feval <= OPTIONS.MaxFEval
    ReplaceFlag = false;
    for q = 1 : OPTIONS.numVar
        Children(q).chrom = Population(1).chrom;
        features = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand;
        Children(q).chrom(q) = features(q);
        Children(q) = OPTIONS.CostFunction(Children(q), OPTIONS);
        if Children(q).cost < Population(1).cost
            Population(1) = Children(q);
            ReplaceFlag = true;
        end        
    end
    feval = feval + OPTIONS.numVar;
    if ~ReplaceFlag
        Population = OPTIONS.InitFunction(OPTIONS);
        Population = OPTIONS.CostFunction(Population, OPTIONS);
        feval = feval + 1;
    end
    if Population(1).cost < Best.cost
        Best = Population(1);
    end
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
    GenIndex = GenIndex + 1;
    fevalArr(GenIndex) = feval;
end
indices = find(fevalArr <= OPTIONS.MaxFEval);
MinCost = MinCost(1 : indices(end));
fevalArr = fevalArr(1 : indices(end));
Conclude(DisplayFlag, OPTIONS, Best, MinCost, AvgCost, MinConstrViol, AvgConstrViol, fevalArr);
return