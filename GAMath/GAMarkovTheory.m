function GAMarkovTheory(pm, pc, n, N, f)

% Using the example on p. 121 of Reeves & Rowe, compute the state transition
% matrix of a simple GA using selection and mutation.
% INPUTS: pm = mutation probability per bit
%         pc = crossover probability
%         n = cardinality of search set
%         N = population size
%         f = n-element vector of fitnesses

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
q = log2(n); % number of bits per chromosome
Chroms = zeros(N, q); % enumerate all possible chromosomes in search space
for i = 1 : n
    temp = i - 1;
    for j = q-1 : -1 : 0
        if temp >= 2^j
            Chroms(i,q-j) = 1;
            temp = temp - 2^j;
        end
    end
end
[~, BestIndex] = max(f); % find the indices of the solutions that are most fit
NumPossible = nchoosek(n+N-1, N); % total number of possible populations
Pops = EnumPops(n, N); % enumerate all possible populations
% Mutation matrix - M(i,j) = M(j,i) = the probability that individual i mutates to individual j
M = zeros(n, n);
for i = 1 : n
    for j = 1 : n
        numbitschange = sum(abs(Chroms(i,:) - Chroms(j,:)));
        M(i,j) = pm^numbitschange * (1-pm)^(q-numbitschange);
        M(j,i) = M(i,j);
    end
end
% Crossover matrix - r(i,j,k) is the probability that x(i) and x(j) cross to form x(k)
r = zeros(n,n,n);
for i = 1 : n
    for j = 1 : n
        offspringcounter = zeros(1, n);
        if i == j
            r(i,i,i) = 1;
        else
            for k = 1 : n
                for m = 1 : q-1
                    offspring1 = [Chroms(i,1:m), Chroms(j, (m+1):q)];
                    offspring2 = [Chroms(j,1:m), Chroms(i, (m+1):q)];
                    if isequal(Chroms(k, :), offspring1)
                        offspringcounter(k) = offspringcounter(k) + 1;
                    end
                    if isequal(Chroms(k, :), offspring2)
                        offspringcounter(k) = offspringcounter(k) + 1;
                    end
                end
                r(i,j,k) = pc * offspringcounter(k) / (2*(q-1));
            end
            r(i,j,i) = r(i,j,i) + 1 - pc;
        end
    end
end
% Compute the transition matrix based on selection and mutation
Q = ones(NumPossible, NumPossible);
disp([num2str(NumPossible), ' population distributions']);
Ps = zeros(1, n);
for j = 1 : NumPossible
    if j > 1, fprintf('\b\b\b\b\b\b\b\b'), end
    fprintf('%5d...', j);
    % For the jth population distribution, compute the probability of selecting individual i
    TotalFitness = sum(Pops(j,:).*f); % Sum of all fitnesses of population distribution j
    for i = 1 : n
        Ps(i) = Pops(j,i) * f(i) / TotalFitness;
    end
    Pcs=zeros(1,n); % probability of obtaining some individual after crossover and selection
    for k = 1:n
        for i = 1:n
            for j1 = 1:n
                Pcs(k) = Pcs(k) + r(i,j1,k) * Ps(i) * Ps(j1);
            end
        end
    end
    Pmcs = zeros(1, n); % probability of obtaining some individual after crossover, selection, mutation
    for p = 1 : n
        Pmcs(p) = sum(Pcs .* M(p,:));
    end
    % For the jth pop distribution, compute the probability of obtaining the kth pop distribution
    for k = 1 : NumPossible
        for i = 1 : n
            if (Pops(k,i) == 0), continue, end;
            Q(j,k) = Q(j,k) * Pmcs(i)^(Pops(k,i)) / factorial(Pops(k,i));
        end
        Q(j,k) = Q(j,k) * factorial(N);
    end
end
fprintf('\n');
% Use Davis-Principe theorem (Theorem 5.4 in Reeves and Rowe) to get the limiting distribution
denom = 0;
for u = 1 : NumPossible
    Qu = Q;
    Qu(:,u) = 0;
    denom = denom + det(Qu - eye(size(Qu)));
end
qvDavis = zeros(1, NumPossible);
for v = 1 : NumPossible
    Qv = Q;
    Qv(:,v) = 0;
    qvDavis(v) = det(Qv - eye(size(Qv))) / denom;
end
% Compute the limiting distribution (not as robust as Davis-Principe)
Qinf = Q^1000;
qvBrute = Qinf(1, :);
% Use the Fundamental Markov Theorem, knowing that Matlab puts the largest eigenvalue first,
% and the largest eigenvalue is equal to 1 (not as robust as Davis-Principe)
[v, ~] = eig(Q');
qvFund = v(:,1) / sum(v(:,1));
disp(['GA: pm = ', num2str(pm), ', pc = ', num2str(pc),', n (search space size) = ', num2str(n), ...
    ', N (pop size) = ', num2str(N)]);
disp(['f = ', num2str(f)]);
% Find the largest probabilities, and the populations corresponding to those probabilities
[yDavis, ndxDavis] = sort(qvDavis, 'descend');
display(['Davis: Prob(', num2str(Pops(ndxDavis(1),:)), ') = ', num2str(yDavis(1))]);
display(['       Prob(', num2str(Pops(ndxDavis(2),:)), ') = ', num2str(yDavis(2))]);
display(['       Prob(', num2str(Pops(ndxDavis(3),:)), ') = ', num2str(yDavis(3))]);
[yFund, ndxFund] = sort(qvFund, 'descend');
display(['Fund:  Prob(', num2str(Pops(ndxFund(1),:)), ') = ', num2str(yFund(1))]);
display(['       Prob(', num2str(Pops(ndxFund(2),:)), ') = ', num2str(yFund(2))]);
display(['       Prob(', num2str(Pops(ndxFund(3),:)), ') = ', num2str(yFund(3))]);
[yBrute, ndxBrute] = sort(qvBrute, 'descend');
display(['Brute: Prob(', num2str(Pops(ndxBrute(1),:)), ') = ', num2str(yBrute(1))]);
display(['       Prob(', num2str(Pops(ndxBrute(2),:)), ') = ', num2str(yBrute(2))]);
display(['       Prob(', num2str(Pops(ndxBrute(3),:)), ') = ', num2str(yBrute(3))]);
% Find the percent of populations that do not have any optimal solutions
NonOptimal = 0;
for i = 1 : NumPossible
    if sum(Pops(i, BestIndex)) == 0
        NonOptimal = NonOptimal + qvDavis(i);
    end
end
display(['Prob(no optimal solutions) = ', num2str(NonOptimal)]);
% Find the percent of uniform optimal populations
AllOptimal = 0;
for i = 1 : NumPossible
    if sum(Pops(i, BestIndex)) == N
        AllOptimal = AllOptimal + qvDavis(i);
    end
end
display(['Prob(uniform optimal population) = ', num2str(AllOptimal)]);