function [MinCost] = SADimension(ProblemFunction, DisplayFlag, GenLimit, betaOdd, betaEven)

% Simulated annealing for minimizing a continuous function.
% Modified to allow different cooling schedules for different dimensions

% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         GenLimit = generation limit
%         betaOdd = cooling parameter for odd-numbered indices: T = T / (1 + beta * T)
%         betaEven = cooling parameter for even-numbered indices: T = T / (1 + beta * T)
% OUTPUT: MinCost = array of best solution, one element for each generation

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @AckleyScaled;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 10000;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('betaOdd', 'var') || isempty(betaOdd)
    betaOdd = 0.001;
end
if ~exist('betaEven', 'var') || isempty(betaEven)
    betaEven = 0.001;
end

% Initialization
OPTIONS.popsize = 1;
OPTIONS.clearDups = false;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);

beta = zeros(1, OPTIONS.numVar);
for i = 1 : 2 : OPTIONS.numVar
    beta(i) = betaOdd;
end
for i = 2 : 2 : OPTIONS.numVar
    beta(i) = betaEven;
end

T0 = 100;
T = T0 * ones(1, OPTIONS.numVar);
BestSoFar = Population(1); % best individual found so far

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    T = T ./ (1 + beta .* T);
    % Generate a candidate solution
    Candidate.chrom = Population(1).chrom + sqrt(T) .* randn(1, OPTIONS.numVar);
    Candidate.chrom = max( min( Candidate.chrom, OPTIONS.MaxDomain ), OPTIONS.MinDomain );
    Candidate = OPTIONS.CostFunction(Candidate, OPTIONS);
    % Decide whether to keep the current solution or replace it with the new candidate solution
    if Candidate.cost < Population(1).cost
        Population(1) = Candidate;
        if Candidate.cost < BestSoFar.cost
            BestSoFar = Candidate;
        end
    else
        Replaced = false;
        for i = 1 : OPTIONS.numVar
            if rand < exp(-1 / T(i))
                Population(1).chrom(i) = Candidate.chrom(i);
                Replaced = true;
            end
        end
        if Replaced
            Population = OPTIONS.CostFunction(Population, OPTIONS);
        end
    end
    DisplayThisTime = DisplayFlag && (mod(GenIndex, 1000) == 0);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayThisTime);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
if DisplayFlag
    xlabel('Iteration')
    ylabel('Best Cost So Far')
end
return
