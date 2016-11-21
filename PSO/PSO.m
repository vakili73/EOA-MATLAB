function [MinCost] = PSO(ProblemFunction, DisplayFlag, c1, c2, c3, KFactor, GenLimit)

% Particle swarm optimization for optimizing a general function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%           the handles of the initialization and cost functions
%         DisplayFlag = whether or not to display information during iterations and plot results
%         c1 = cognitive constant
%         c2 = social constant for neighborhood interaction
%         c3 = social constant for global interaction
%         KFactor = constriction coefficient
%         GenLimit = generation limit

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('c1', 'var') || isempty(c1)
    c1 = 2.1;
end
if ~exist('c2', 'var') || isempty(c2)
    c2 = 2.1;
end
if ~exist('c3', 'var') || isempty(c3)
    c3 = 2.1;
end
if ~exist('KFactor', 'var') || isempty(KFactor)
    KFactor = 0.9;
end
Constriction = KFactor * 2 / (c1 + c2 + c3 - 2);
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 50;
end
OPTIONS.Maxgen = GenLimit;

OPTIONS.numVar = 20;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
OPTIONS.neighbors = 5;
OPTIONS.neighbors = OPTIONS.neighbors + 1; % size of particle swarm neighborhood

for k = 1 : OPTIONS.popsize
    Population(k).vel(1 : OPTIONS.numVar) = 0; % initialize velocities
end
pbest = Population; % personal best of each particle
nbest = Population; % neighborhood best of each particle
gbest = Population(1); % global best

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Update the global best if needed
    if Population(1).cost < gbest.cost
        gbest = Population(1);
    end
    % Update the personal best and neighborhood best for each particle
    for i = 1 : OPTIONS.popsize 
        % Update each personal best as needed
        if Population(i).cost < pbest(i).cost
            pbest(i) = Population(i);
        end
        % Update each neighborhood best as needed
        Distance = zeros(OPTIONS.popsize, 1);
        for j = 1 : OPTIONS.popsize 
            Distance(j) = norm(Population(i).chrom - Population(j).chrom);
        end
        [~, indices] = sort(Distance);
        nbest(i).cost = inf;
        for j = 2 : OPTIONS.neighbors
            nindex = indices(j);
            if Population(nindex).cost < nbest(i).cost
                nbest(i) = Population(nindex);
            end
        end
    end
    % Update the position and velocity of each particle (except the elites)
    for i = OPTIONS.Keep+1 : OPTIONS.popsize
        r = rand(3, OPTIONS.numVar);
        x = Population(i).chrom;
        deltaVpersonal = c1 * r(1,:) .* (pbest(i).chrom - x);
        deltaVneighborhood = c2 * r(2,:) .* (nbest(i).chrom - x);
        deltaVglobal = c3 * r(3,:) .* (gbest.chrom - x);
        Population(i).vel = Constriction * (Population(i).vel + deltaVpersonal + deltaVneighborhood + deltaVglobal);
        Population(i).chrom = x + Population(i).vel;
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