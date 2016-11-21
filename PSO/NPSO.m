function [MinCost] = NPSO(ProblemFunction, DisplayFlag, c1, c2, c3, c4, c5, c6, KFactor, GenLimit)

% Negative reinforcment PSO for optimizing a general function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%           the handles of the initialization and cost functions
%         DisplayFlag = whether or not to display information during iterations and plot results
%         c1 = cognitive constant
%         c2 = social constant for neighborhood interaction
%         c3 = social constant for global interaction
%         c4 = negative reinforcement cognitive constant
%         c5 = negative reinforcement social constant for neighborhood interaction
%         c3 = negative reinforcement social constant for global interaction
%         KFactor = constriction coefficient
%         GenLimit = generation limit
%         PopSize = population size

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('c1', 'var')
    c1 = 2.1;
end
if ~exist('c2', 'var')
    c2 = 2.1;
end
if ~exist('c3', 'var')
    c3 = 2.1;
end
if ~exist('c4', 'var')
    c4 = 0.0;
end
if ~exist('c5', 'var')
    c5 = 0.0;
end
if ~exist('c6', 'var')
    c6 = 0.0;
end
if ~exist('KFactor', 'var')
    KFactor = 0.9;
end
Constriction = KFactor * 2 / (c1 + c2 + c3 - 2);
if ~exist('GenLimit', 'var')
    GenLimit = 50;
end
OPTIONS.Maxgen = GenLimit;

OPTIONS.popsize = 20;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
OPTIONS.neighbors = 5; % size of particle swarm neighborhood

for k = 1 : OPTIONS.popsize
    Population(k).vel(1 : OPTIONS.numVar) = 0; % initialize velocities
end
pbest = Population; % personal best of each particle
pworst = Population; % personal worst of each particle
nbest = Population; % neighborhood best of each particle
nworst = Population; % neighborhood worst of each particle
gbest = Population(1); % global best
gworst = Population(end); % global worst

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Update the global best and worst if needed
    if Population(1).cost < gbest.cost
        gbest = Population(1);
    elseif Population(end).cost > gworst.cost
        gworst = Population(end);
    end
    % Update the personal best and neighborhood best for each particle
    for i = 1 : OPTIONS.popsize 
        % Update each personal best and worst as needed
        if Population(i).cost < pbest(i).cost
            pbest(i) = Population(i);
        elseif Population(i).cost > pworst(i).cost
            pworst(i) = Population(i);
        end
        % Update each neighborhood best and worst as needed
        Distance = zeros(OPTIONS.popsize, 1);
        for j = 1 : OPTIONS.popsize 
            Distance(j) = norm(Population(i).chrom - Population(j).chrom);
        end
        [~, indices] = sort(Distance);
        nbest(i).cost = inf;
        nworst(i).cost = -inf;
        for j = 2 : OPTIONS.neighbors
            nindex = indices(j);
            if Population(nindex).cost < nbest(i).cost
                nbest(i) = Population(nindex);
            elseif Population(nindex).cost > nworst(i).cost
                nworst(i) = Population(nindex);
            end
        end
    end
    % Update the position and velocity of each particle (except the elites)
    for i = OPTIONS.Keep+1 : OPTIONS.popsize
        r = rand(6, OPTIONS.numVar);
        x = Population(i).chrom;
        deltaVpersonal = c1 * r(1,:) .* (pbest(i).chrom - x);
        deltaVneighborhood = c2 * r(2,:) .* (nbest(i).chrom - x);
        deltaVglobal = c3 * r(3,:) .* (gbest.chrom - x);
        deltaVpworst = -c4 * r(4,:) .* (pworst(i).chrom - x);
        deltaVnworst = -c5 * r(5,:) .* (nworst(i).chrom - x);
        deltaVgworst = -c6 * r(6,:) .* (gworst.chrom - x);
        Population(i).vel = Constriction * (Population(i).vel + deltaVpersonal + deltaVneighborhood + ...
             + deltaVglobal + deltaVpworst + deltaVnworst + deltaVgworst);
        Population(i).chrom = x + Population(i).vel;
    end 
    % Limit each solution to the search domain
    for i = 1 : OPTIONS.popsize
        Population(i).chrom = min(max([Population(i).chrom], OPTIONS.MinDomain), OPTIONS.MaxDomain);
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