function GADyn2

% Dynamic system model of simple selection-mutation GA.

f = [5 2 2 3 2 3 3 4]; % fitnesses
mu = 0.02; % mutation rate
genLimit = 100; % number of generations to run
PopSize = 500; % population size
Chroms = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
pInit = [0 0 0 .1 0 .1 .1 .7]';
% Mutation matrix - M(j,i)=M(i,j) is the probability that island j mutates to island i
M = zeros(length(f), length(f));
for i = 1 : length(f)
    for j = 1 : length(f)
        d = sum(abs(Chroms(i,:) - Chroms(j,:))); % Hamming distance
        M(i,j) = mu^d * (1 - mu)^(size(Chroms,2) - d);
    end
end
[v, ~] = eig(M * diag(f));
for i = 1 : length(f)
    v(:,i) = v(:,i) / sum(v(:,i));
    disp(['eigenvector # ', num2str(i), ' = ', num2str(v(:,i)')]);
end
% Create the initial population
x = zeros(PopSize, size(Chroms,2));
LastEnd = 0;
for i = 1 : length(pInit)
    x(LastEnd+1:LastEnd+pInit(i)*PopSize, :) = ones(pInit(i)*PopSize, 1) * Chroms(i, :);
    LastEnd = LastEnd + pInit(i) * PopSize;
end
% Initialize arrays
pArrSim = zeros(genLimit+1, length(pInit));
pArrSim(1, :) = pInit * PopSize;
HammingMetaArr = zeros(genLimit, 1);
HammingStableArr = zeros(genLimit, 1);
fitness = zeros(PopSize, 1);
selected = zeros(PopSize, size(Chroms, 2));
for gen = 1 : genLimit
    % Calculate fitness
    for pop = 1 : PopSize
        for i = 1 : length(f)
            if isequal(x(pop,:), Chroms(i,:))
                fitness(pop) = f(i);
                break
            end
        end
    end
    % Selection
    for pop = 1 : PopSize
        partsum = 0;
        j = 0;
        randomnumber = rand * sum(fitness);
        while (partsum <= randomnumber) && (j <= PopSize)
            j = j + 1;
            partsum = partsum + fitness(j);
        end
        selected(pop, :) = x(j, :);
    end
    x = selected;
    % Mutation
    for pop = 1 : PopSize
        for i = 1 : size(x, 2)
            if rand < mu
                x(pop, i) = ~x(pop, i);
            end
        end
    end
    % Compute distance of population from metastable and stable fixed points
    HammingMeta = 0;
    HammingStable = 0;
    for pop = 1 : PopSize
        HammingMeta = HammingMeta + sum(abs(x(pop,:) - Chroms(1,:)));
        HammingStable = HammingStable + sum(abs(x(pop,:) - Chroms(8,:)));
        for i = 1 : length(pInit)
            if isequal(x(pop,:), Chroms(i,:))
                pArrSim(gen+1, i) = pArrSim(gen+1, i) + 1;
            end
        end
    end
    HammingMetaArr(gen) = HammingMeta;
    HammingStableArr(gen) = HammingStable;
end
close all; 
figure; set(gca, 'fontsize', 14); hold on; box on;
gen = 1 : genLimit;
plot(gen, HammingMetaArr/PopSize, 'b:', 'linewidth', 2);
plot(gen, HammingStableArr/PopSize, 'r-', 'linewidth', 2);
legend('Distance from [0 0 0]', 'Distance from [1 1 1]');
xlabel('Generation'); ylabel('Hamming distance');

figure; set(gca, 'fontsize', 14); hold on; box on;
gen = 0 : genLimit;
plot(gen, pArrSim(:, 1)/PopSize, 'b:', 'linewidth', 2);
plot(gen, pArrSim(:, 8)/PopSize, 'r-', 'linewidth', 2);
legend('p_1', 'p_8');
xlabel('Generation'); ylabel('Proportion');

% Analytical results
% Lambda = d;
% d = diag(Lambda);
% q = inv(v) * pInit;
% pArr = zeros(length(pInit), genLimit+1);
% pArr(:, 1) = pInit;
% for gen = 1 : genLimit
%     q = Lambda * q / (d' * q);
%     pArr(:, gen+1) = v * q;
% end

pArr(:, 1) = pInit;
p = pInit;
for gen = 1 : genLimit
    p = M' * diag(f) * p / (f * p);
    pArr(:, gen+1) = p;
end


figure; set(gca, 'fontsize', 14); hold on; box on;
gen = 0 : genLimit;
plot(gen, pArr(1,:), 'b:', 'linewidth', 2);
plot(gen, pArr(8,:), 'r-', 'linewidth', 2);
legend('Proportion of [0 0 0]', 'Proportion of [1 1 1]');
xlabel('Generation'); ylabel('Proportion');