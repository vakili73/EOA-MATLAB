function [MinCost, AvgCost] = GA(ProblemFunction, DisplayFlag, RandSeed, GenLimit, GrayFlag, NumElites, ...
    PopSize, Pmutate, MinDomain, MaxDomain, Dimension, ElitismOption)
% Genetic algorithm for optimizing a general function.
% INPUTS: ProblemFunction = the handle of the function that returns 
%            the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = whether or not to display information during iterations and plot results.
%         RandSeed = random number seed
%         GenLimit = number of generations
%         GrayFlag = whether to use gray coding or binary coding
%         NumElites = number of elite individuals each generation
%         PopSize = population size
%         Pmutate = mutation rate per bit (discrete problems) or per solution feature (continuous problems)
%         MinDomain = lower bound of search domain
%         MaxDomain = upper bound of search domain
%         Dimension = problem dimension
%         ElitismOption = 1 means generate (N-E) children, where N = population size, and E = number of elites
%                         2 means generate N children and replace the worst E with the elites
% OUTPUTS: MinCost = array of best cost values, one element per generation
%          AvgCost = array of average cost values, one element per generation
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @AckleyDisc;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = round(sum(100*clock));
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 10;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('GrayFlag', 'var') || isempty(GrayFlag)
    GrayFlag = false;
end
OPTIONS.Gray = GrayFlag;
if ~exist('NumElites', 'var') || isempty(NumElites)
    NumElites = 2;
end
OPTIONS.Keep = NumElites;
if ~exist('PopSize', 'var') || isempty(PopSize)
    PopSize = 20;
end
OPTIONS.popsize = PopSize;
if ~exist('Pmutate', 'var') || isempty(Pmutate)
    Pmutate = 0.02;
end
OPTIONS.pmutate = Pmutate;
if ~exist('Dimension', 'var') || isempty(Dimension)
    Dimension = 2;
end
OPTIONS.numVar = Dimension;
if ~exist('MinDomain', 'var') || isempty(MinDomain)
    MinDomain = -5;
end
OPTIONS.MinDomain = MinDomain * ones(1, Dimension);
if ~exist('MaxDomain', 'var') || isempty(MaxDomain)
    MaxDomain = 5;
end
OPTIONS.MaxDomain = MaxDomain * ones(1, Dimension);
if ~exist('ElitismOption', 'var') || isempty(ElitismOption)
    ElitismOption = 1;
end
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);
Xover_Type = 1; % crossover type: 1 = single point, 2 = two point, 3 = uniform
pcross = 1; % crossover probability
if ElitismOption == 1
    child = zeros(OPTIONS.popsize-OPTIONS.Keep, OPTIONS.numVar);
else
    child = zeros(OPTIONS.popsize, OPTIONS.numVar);
end
if DisplayFlag && isfield(OPTIONS, 'Dim') && (OPTIONS.Dim == 2)
    PlotContour(OPTIONS, Population, 1);
end
% Begin the evolution loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Save the elites
    Elites = Population(1 : OPTIONS.Keep);
    % Compute fitness, which increases as cost decreases.
    Fitness = -[Population.cost];
    Fitness = Fitness - min(Fitness);
    Fitness = Fitness + max(Fitness);
    FitnessSum = sum(Fitness);
    for k = 1 : 2 : size(child, 1) % begin selection/crossover loop
        % Select two parents to mate and create two children - roulette wheel selection
        mate = [0 0];
        for selParents = 1 : 2
            Random_Cost = rand * FitnessSum;
            Select_Cost = Fitness(1);
            Select_index = 1;
            while Select_Cost < Random_Cost 
                Select_index = Select_index + 1;
                if Select_index >= OPTIONS.popsize, break, end
                Select_Cost = Select_Cost + Fitness(Select_index);
            end
            mate(selParents) = Select_index;
        end
        Parent(1, :) = Population(mate(1)).chrom;
        Parent(2, :) = Population(mate(2)).chrom;
        % Crossover
        switch Xover_Type
            case 1
                % single-point crossover
                if rand < pcross
                    Xover_Pt = ceil(rand * OPTIONS.numVar);
                    child(k, :) = [Parent(1, 1:Xover_Pt), Parent(2, Xover_Pt+1:OPTIONS.numVar)];
                    child(k+1, :) = [Parent(2, 1:Xover_Pt), Parent(1, Xover_Pt+1:OPTIONS.numVar)];
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
            case 2
                % two-point crossover
                if rand < pcross
                    Xover_Pt1 = ceil(rand * OPTIONS.numVar);
                    Xover_Pt2 = ceil(rand * OPTIONS.numVar);
                    if Xover_Pt1 > Xover_Pt2
                        temp = Xover_Pt2;
                        Xover_Pt2 = Xover_Pt1;
                        Xover_Pt1 = temp;
                    end
                    child(k, :) = [Parent(1, 1:Xover_Pt1) Parent(2, Xover_Pt1+1:Xover_Pt2) Parent(1, Xover_Pt2+1:OPTIONS.numVar)];
                    child(k+1, :) = [Parent(2, 1:Xover_Pt1) Parent(1, Xover_Pt1+1:Xover_Pt2) Parent(2, Xover_Pt2+1:OPTIONS.numVar)];
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
            case 3 
                % uniform crossover
                if rand < pcross
                    for i = 1 : OPTIONS.numVar
                        if rand < 0.5
                            child(k, i) = Parent(1, i);
                            child(k+1, i) = Parent(2, i);
                        else
                            child(k, i) = Parent(2, i);
                            child(k+1, i) = Parent(1, i);
                        end
                    end
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
        end
    end % end selection/crossover loop
    % Mutation - but don't allow the elites to be mutated
    for ind = 1 : size(child, 1)
        for parnum = 1 : OPTIONS.numVar
            if rand < OPTIONS.pmutate
                if ~isempty(strfind(func2str(ProblemFunction), 'Disc'))
                    % The problem is binary, so perform binary mutation
                    child(ind, parnum) = ~child(ind, parnum);
                else
                    % The problem is continuous, so perform uniform continuous mutation
                    child(ind, parnum) = OPTIONS.MinDomain(parnum) + ...
                        rand * (OPTIONS.MaxDomain(parnum) - OPTIONS.MinDomain(parnum));
                end
            end
        end
    end
    % Create a new population for the next generation comprised of this generation's children combined
    % with the previous generation's elites
    for ind = 1 : size(child, 1)
        Population(ind).chrom = child(ind, :);
    end
    for ind = 1 : OPTIONS.Keep
        Population(size(child,1)+ind).chrom = Elites(ind).chrom;
    end
    % Make sure the population does not have duplicates.
    if OPTIONS.clearDups
        Population = ClearDups(Population, OPTIONS);
    end
    if isfield(OPTIONS, 'FeasibleFunction') && ~isempty(OPTIONS.FeasibleFunction)
        % Make sure each individual is legal.
        Population = OPTIONS.FeasibleFunction(Population, OPTIONS);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    Population = PopSort(Population, OPTIONS.ConstrMethod, GenIndex);
    Population = Population(1 : OPTIONS.popsize);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
    if DisplayFlag && isfield(OPTIONS, 'Dim') && (OPTIONS.Dim == 2)
        if (GenIndex == 4)
            PlotContour(OPTIONS, Population, 2);
        elseif (GenIndex == 7)
            PlotContour(OPTIONS, Population, 3);
        elseif (GenIndex == 10)
            PlotContour(OPTIONS, Population, 4);
        end
    end
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return