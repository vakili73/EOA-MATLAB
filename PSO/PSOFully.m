function [MinCost] = PSOFully(ProblemFunction, DisplayFlag, phiMax, KFactor, GenLimit)

% Fully informed particle swarm optimization for optimizing a general function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%            the handles of the initialization and cost functions
%         DisplayFlag = whether or not to display information during iterations and plot results
%         phiMax = maximum value of each component of phi
%         KFactor = constriction coefficient3
%         GenLimit = generation limit

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('phiMax', 'var') || isempty(phiMax)
    phiMax = 2;
end
if ~exist('KFactor', 'var') || isempty(KFactor)
    KFactor = 0.9;
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 40;
end
OPTIONS.Maxgen = GenLimit;
OPTIONS.popsize = 50;

OPTIONS.numVar = 20;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
Constriction = KFactor * 2 / (3 * phiMax - 2);

for k = 1 : OPTIONS.popsize
    Population(k).vel(1 : OPTIONS.numVar) = 0; % initialize velocities
end
BestInd = Population; % personal best of each particle
Distance = zeros(OPTIONS.popsize, OPTIONS.popsize); % distance between particles

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Update the personal best for each particle
    for i = 1 : OPTIONS.popsize
        if Population(i).cost < BestInd(i).cost
            BestInd(i) = Population(i);
        end
        % Update the distance array
        for j = 1 : OPTIONS.popsize 
            Distance(i, j) = norm(Population(i).chrom - Population(j).chrom);
        end
    end
    % Update the position and velocity of each particle (except the elites)
    WorstCost = max([Population.cost]);
    for i = OPTIONS.Keep+1 : OPTIONS.popsize
        MaxDistance = max(Distance(i, :));
        CostFactors = WorstCost - [Population.cost];
        DistanceFactors = MaxDistance - Distance(i, :);
        
        % ScaleFactor = max(CostFactors) / max(DistanceFactors); This was corrected on Aug. 4, 2014, to be consistent with Equation (11.35)
        ScaleFactor = max(CostFactors) / MaxDistance;
        
        Weights = CostFactors + ScaleFactor * DistanceFactors;
        %Weights = Weights / sum(Weights);
        phi = phiMax * rand(1, OPTIONS.popsize);
        x = reshape([BestInd.chrom], OPTIONS.popsize, OPTIONS.numVar);
        psubi = (Weights .* phi) * x / (Weights * phi');
        Population(i).vel = Constriction * (Population(i).vel + mean(phi) * (psubi - Population(i).chrom));
        Population(i).chrom = Population(i).chrom + Population(i).vel;
        Population(i).chrom = max( min(Population(i).chrom, OPTIONS.MaxDomain), OPTIONS.MinDomain);
    end 
    if OPTIONS.clearDups
        Population = ClearDups(Population, OPTIONS);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    Population = PopSort(Population);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return