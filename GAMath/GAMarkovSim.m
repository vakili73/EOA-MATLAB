function [OutputPops, OutputProbs] = GAMarkovSim(pm, pc, n, N, f, genLimit)

% Using the example on p. 121 of Reeves & Rowe, run a simple GA and 
% plot the proportion of various populations
% INPUT: pm = mutation probability per bit
%        pc = crossover probability
%        n = cardinality of search set
%        N = population size
%        f = n-element vector of fitnesses
%        genLimit = number of generations to run
% OUTPUTS: OutputPops = 10 most likely populations
%          OutputProbs = 10 highest population probabilities

if ~exist('pm', 'var') || isempty(pm)
    pm = 0.1;
end
if ~exist('pc', 'var') || isempty(pc)
    pc = 0.9;
end
if ~exist('n', 'var') || isempty(n)
    if exist('f', 'var') && ~isempty(f)
        n = length(f);
    else
        n = 8;
    end
end
if ~exist('N', 'var') || isempty(N)
    N = 3;
end
if ~exist('f', 'var')
    f = [1 2 2 3 2 3 3 4]; % one-max problem (Example 4.9)
    f = [5 2 2 3 2 3 3 4]; % deceptive problem (Example 4.10)
end
if ~exist('genLimit', 'var')
    genLimit = 20000;
end
numbits = log2(n); % number of bits per chromosome
Chroms = zeros(N, numbits); % enumerate all possible chromosomes in search space
for i = 1 : n
    temp = i - 1;
    for j = numbits-1 : -1 : 0
        if temp >= 2^j
            Chroms(i,numbits-j) = 1;
            temp = temp - 2^j;
        end
    end
end
[~, BestIndex] = max(f); % find the indices of the solutions that are most fit
x = round(rand(N, numbits)); % initial population
NumPossible = nchoosek(n+N-1, N); % total number of possible populations
PopCount = zeros(1, NumPossible); % Running count of each population possibility
Pct = zeros(genLimit, NumPossible); % Running percent of each population distribution
NoOptimal = 0; % Running count of number of populations containing no optimal solutions
NoOptimalPct = zeros(1, genLimit); % Running percentage of number of populations containing no optimal solutions
AllOptimal = 0; % Running count of number of uniform optimal populations
AllOptimalPct = zeros(1, genLimit); % Running percentage of number of uniform optimal populations
Pops = EnumPops(n, N); % enumerate all possible populations
disp([num2str(genLimit), ' generations']);
fitness = zeros(1, N);
selected = zeros(2*N, numbits);
x1 = zeros(size(x));
for gen = 1 : genLimit
    if mod(gen,1000) == 0
        if gen > 1000, fprintf('\b\b\b\b\b\b\b\b\b\b'), end
        fprintf('%7d...', gen);
    end
    % Calculate fitness
    for pop = 1 : N
        for i = 1 : n
            if isequal(x(pop,:), Chroms(i,:))
                fitness(pop) = f(i);
                break
            end
        end
    end
    % Selection
    for pop = 1 : 2*N
        partsum = 0;
        j = 0;
        randomnumber = rand * sum(fitness);
        while (partsum <= randomnumber) && (j <= N)
            j = j + 1;
            partsum = partsum + fitness(j);
        end
        selected(pop, :) = x(j, :);
    end
    x = selected;
    % Mutation
    for pop = 1 : 2*N
        for i = 1 : size(x, 2)
            if rand < pm
                x(pop, i) = ~x(pop, i);
            end
        end
    end
    % Crossover
    for pop = 1 : N
        parent1 = pop;
        if rand < pc
            parent2 = parent1 + N;
            crossoverpoint = min(numbits-1, floor((numbits-1)*rand) + 1);
            if rand < 0.5
                x1(pop,:) = [x(parent1,1:crossoverpoint), x(parent2,crossoverpoint+1:numbits)];
            else
                x1(pop,:) =  [x(parent2,1:crossoverpoint), x(parent1,crossoverpoint+1:numbits)];
            end
        else
            x1(pop,:) = x(parent1, :);
        end
    end
    x = x1;
    % Count populations containing each possible distribution of individuals
    ChromCount = zeros(1, n);
    for i = 1 : N
        for j = 1 : n
            if isequal(x(i,:), Chroms(j,:))
                ChromCount(j) = ChromCount(j) + 1;
                break
            end
        end        
        for k = 1 : NumPossible
            if isequal(ChromCount, Pops(k, :))
                PopCount(k) = PopCount(k) + 1;
                break
            end
        end
    end
    Pct(gen, :) = PopCount / gen;
    % Find the running percentage of populations that do not have any optimal solutions,
    % the the running percentage of populations that are uniformly optimal.
    NumOptimaThisGen = sum(ChromCount(BestIndex));
    if NumOptimaThisGen == 0
        NoOptimal = NoOptimal + 1;
    else
        if NumOptimaThisGen == N
            AllOptimal = AllOptimal + 1;
        end
    end
    NoOptimalPct(gen) = NoOptimal / gen;
    AllOptimalPct(gen) = AllOptimal / gen;
end
fprintf('\n');
disp(['GA: pm = ', num2str(pm), ', pc = ', num2str(pc),', n (search space size) = ', num2str(n), ...
    ', N (pop size) = ', num2str(N)]);
disp(['f = ', num2str(f)]);
% Display and plot the most likely population distributions
[y, ndx] = sort(PopCount, 'descend');
display(['Prob(', num2str(Pops(ndx(1),:)), ') = ', num2str(y(1) / genLimit)]);
display(['Prob(', num2str(Pops(ndx(2),:)), ') = ', num2str(y(2) / genLimit)]);
display(['Prob(', num2str(Pops(ndx(3),:)), ') = ', num2str(y(3) / genLimit)]);
display(['Prob(no optimal solutions) = ', num2str(NoOptimalPct(end))]);
display(['Prob(all optimal solutions) = ', num2str(AllOptimalPct(end))]);
OutputPops = Pops(ndx(1:10),:);
OutputProbs = [y(1:10) / genLimit, NoOptimalPct(end), AllOptimalPct(end)];
close all; figure; set(gca, 'fontsize', 14); hold on; box on;
plot(100*NoOptimalPct, 'm:', 'linewidth', 2);
plot(100*Pct(:,ndx(1)), 'b-.', 'linewidth', 2);
plot(100*Pct(:,ndx(2)), 'r--', 'linewidth', 2); 
plot(100*Pct(:,ndx(3)), 'k-', 'linewidth', 2);
legend('No optimal individuals', num2str(Pops(ndx(1),:)), num2str(Pops(ndx(2),:)), num2str(Pops(ndx(3),:)));
xlabel('generation number');
ylabel('cumlative percent of population');
%title(['Genetic Algorithm, pm = ', num2str(pm)]);