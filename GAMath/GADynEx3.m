function [ProbTheory, ProbSim] = GADynEx3(pm, pc, PopSize, fitness, ProbInit, genLimit, DisplayFlag)

% Dynamic systems analysis of simple GA with single-point crossover (GA/sp).
% Reeves & Rowe, p. 153.

close all
if ~exist('pm', 'var') || isempty(pm)
    pm = 0.01; % mutation probability
end
if ~exist('pc', 'var') || isempty(pc)
    pc = 0.9; % probability of crossover
end
if ~exist('PopSize', 'var') || isempty(PopSize)
    PopSize = 1000; % population size
end
if ~exist('fitness', 'var') || isempty(fitness)
    fitness = [1 2 2 3 2 3 3 4]; % fitness values
end
SearchCardinality = length(fitness); % search space cardinality
if ~exist('ProbInit', 'var') || isempty(ProbInit)
    ProbInit = [0.8; 0.1; 0.1; zeros(SearchCardinality-3, 1)]; % initial population proportions
    %ProbInit = [0.0; 0.1; 0.1; zeros(SearchCardinality-4, 1); 0.8]; % initial population proportions
end
if ~exist('genLimit', 'var') || isempty(genLimit)
    genLimit = 80; % number of generations to run
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true; % whether or not to plot results
end
[~, BestIndex] = max(fitness); % find the indices of the solutions that are most fit
[~, WorstIndex] = min(fitness); % find the indices of the solutions that are least fit
q = round(log2(SearchCardinality)); % number of bits per individual
Chroms = zeros(SearchCardinality, q); % enumerate all possible solutions (chromosomes)
for i = 1 : SearchCardinality
    temp = i - 1;
    for j = q-1 : -1 : 0
        if temp >= 2^j
            Chroms(i, q-j) = 1;
            temp = temp - 2^j;
        end
    end
end
% Mutation matrix - M(i,j) = M(j,i) is the probability that island j mutates to island i
M = zeros(length(fitness));
for i = 1 : length(fitness)
    for j = 1 : length(fitness)
        HamDist = sum(abs(Chroms(i,:) - Chroms(j,:))); % Hamming distance
        M(i,j) = pm^HamDist * (1 - pm)^(size(Chroms,2) - HamDist);
    end
end
% Initialize the population
v = PopSize * ProbInit; % initial population count
ndx = 0;
x = zeros(PopSize, q);
for i = 1 : SearchCardinality
    for j = ndx+1 : ndx+v(i)
        x(j, :) = Chroms(i, :);
    end
    ndx = ndx + v(i);
end
% Initialize for GA simulation
n = SearchCardinality;
ProbArr = ProbInit;
PopFitness = zeros(1, PopSize);
selected = zeros(2*PopSize, q);
x = zeros(2*PopSize, q);
% Simulate the GA
for gen = 1 : genLimit
    % Calculate fitness
    for pop = 1 : PopSize
        for i = 1 : length(fitness)
            if isequal(x(pop,:), Chroms(i,:))
                PopFitness(pop) = fitness(i);
                break
            end
        end
    end
    % Selection
    for pop = 1 : 2*PopSize
        partsum = 0;
        j = 0;
        randomnumber = rand * sum(PopFitness);
        while (partsum <= randomnumber)
            j = j + 1;
            partsum = partsum + PopFitness(j);
        end
        j = min(j, PopSize);
        selected(pop, :) = x(j, :);
    end
    % Crossover
    for pop = 1 : 2*PopSize
        if rand < pc
            parentsNdx = ceil(2 * PopSize * rand(2, 1));
            parentsNdx = max(parentsNdx, 1); % just in case the random number generator returned 0
            parent1 = selected(parentsNdx(1), :);
            parent2 = selected(parentsNdx(2), :);
            el = ceil((q - 1) * rand); % crossover point
            el = max(el, 1);
            x(pop, :) = [parent1(1:el), parent2(el+1 : q)];
        else
            childNdx = ceil(2 * PopSize * rand);
            childNdx = max(childNdx, 1); % just in case the random number generator returned 0
            x(pop, :) = selected(childNdx, :);
        end
    end
    % Mutation
    for pop = 1 : PopSize
        for i = 1 : size(x, 2)
            if rand < pm
                x(pop, i) = ~x(pop, i);
            end
        end
    end
    % Count the number of each individual
    ChromCount = zeros(1, n);
    for i = 1 : PopSize
        for j = 1 : n
            if isequal(x(i,:), Chroms(j,:))
                ChromCount(j) = ChromCount(j) + 1;
                break
            end
        end
    end
    ProbArr(:, gen+1) = ChromCount / PopSize;
end
ProbSim = mean(ProbArr(:, 10:end), 2);
% Analytical results
% Crossover matrix
r = zeros(n, n, n);
for i = 1 : n
    for j1 = 1 : n
        if i ~= j1 % if an individual is not crossing over with itself
            offspringcounter = zeros(1, n);
            for el = 1 : q - 1
                offspring1 = [Chroms(i,1:el), Chroms(j1,(el+1):q)];
                offspring2 = [Chroms(j1,1:el), Chroms(i,(el+1):q)];
                for k = 1 : n
                    if isequal(offspring1, Chroms(k, :))
                        offspringcounter(k) = offspringcounter(k) + 1;
                    end
                    if isequal(offspring2, Chroms(k, :))
                        offspringcounter(k) = offspringcounter(k) + 1;
                    end
                end
            end
            r(i, j1, :) = pc * offspringcounter / sum(offspringcounter);
            r(i, j1, i) = r(i, j1, i) + (1 - pc) / 2;
            r(i, j1, j1) = r(i, j1, j1) + (1 - pc) / 2;
        else
            r(i,j1,i) = 1;
        end
    end
end
Prob = ProbInit;
ProbTheoryArr = Prob;
ProbElements = zeros(n, 1);
for gen = 1 : genLimit
    for k = 1 : n
        R = squeeze(r(:, :, k));
        ProbElements(k) = Prob' * diag(fitness) * M * R * M' * diag(fitness) * Prob / (fitness * Prob)^2;
    end
    Prob = ProbElements;
    ProbTheoryArr(:, gen+1) = Prob;
end
ProbTheory = Prob;
if DisplayFlag
    gen = 0 : genLimit;
    % Plot percent of best individuals
    figure; hold on;
    set(gca, 'FontSize', 14); set(gca, 'Box', 'on'); set(gca, 'DefaultLineLineWidth', 2);
    for i = 1 : length(BestIndex)
        plot(gen, 100*ProbArr(BestIndex(i), :), 'b-')
        plot(gen, 100*ProbTheoryArr(BestIndex(i), :), 'r--')        
    end
    legend('simulation', 'theory')
    xlabel('generation'); ylabel('percent of optimum');
    title(['GA with single point crossover - pm = ', num2str(pm)]);
    % Plot percent of worst individuals
    figure; hold on;
    set(gca, 'FontSize', 14); set(gca, 'Box', 'on'); set(gca, 'DefaultLineLineWidth', 2);
    for i = 1 : length(WorstIndex)
        plot(gen, 100*ProbArr(WorstIndex(i), :), 'b-')
        plot(gen, 100*ProbTheoryArr(WorstIndex(i), :), 'r--')        
    end
    legend('simulation', 'theory')
    xlabel('generation'); ylabel('percent of least fit');
end