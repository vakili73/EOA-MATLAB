function [MinCost] = SA(ProblemFunction, DisplayFlag, Random, GenLimit, RestartCount, alpha, beta, T0, cross, deriv)

% Simulated annealing for minimizing a continuous function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         GenLimit = generation limit
%         RestartCount = how often to restart annealing if there is no improvement
%         alpha = cooling parameter: T = alpha * T
%         beta = cooling parameter: T = T / (1 + beta * T)
%         T0 = initial temperature
% OUTPUT: MinCost = array of best solution, one element for each generation

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('Random', 'var') || isempty(DisplayFlag)
    Random = 'randn()';
%     Random = 'rand()-rand()';
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 500;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('RestartCount', 'var') || isempty(RestartCount)
%     RestartCount = inf;
    RestartCount = 5;
end
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 0;
end
if ~exist('beta', 'var') || isempty(beta)
    beta = 0.001;
end
if ~exist('T0', 'var') || isempty(T0)
    T0 = 100; 
end
if ~exist('cross', 'var') || isempty(cross)
    cross = 0.05; 
end
% if ~exist('deriv', 'var') || isempty(deriv)
%     deriv = 5;
% end

% Initialization
OPTIONS.popsize = 50;
OPTIONS.clearDups = false;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);

T = T0;
% TGama = T0;
TArray = zeros(OPTIONS.Maxgen+1, 1);
TArray(1) = T;
ConsecFails = 0; % number of generations with no improvement
% BestSoFar = Population(1); % best individual found so far

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    if alpha > 0
        T = alpha * T;
    else
        T = T / (1 + beta * T);
    end
    % Generate a candidate solution
    Candidate = struct();
    for i = 1 : length(Population)
        Candidate(i).chrom = Population(i).chrom + sqrt(T) * eval(RandomInterpreter(Random, '1, OPTIONS.numVar'));
%         Candidate(i).chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* eval(RandomInterpreter(Random, '1, OPTIONS.numVar'));
        Candidate(i).chrom = max( min( Candidate(i).chrom, OPTIONS.MaxDomain ), OPTIONS.MinDomain );
    end
%     Candidate.chrom = Population(1).chrom + sqrt(T) * eval(RandomInterpreter(Random, '1, OPTIONS.numVar'));
%     Candidate.chrom = Population(1).chrom + sqrt(T) * (Random(1, OPTIONS.numVar) - Random(1, OPTIONS.numVar));
%     Candidate.chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
%     Candidate.chrom = max( min( Candidate.chrom, OPTIONS.MaxDomain ), OPTIONS.MinDomain );
    Candidate = OPTIONS.CostFunction(Candidate, OPTIONS);
    % Decide whether to keep the current solution or replace it with the new candidate solution
    Candidate = SortByCost(Candidate);
    if Candidate(1).cost < Population(1).cost
        Population = PopAppend(Population, Candidate);
        Population = SortByCost(Population);
        Population = PopPick(Population, OPTIONS.popsize);
%         ConsecFails = 0;
%         if BestCandidate.cost < BestSoFar.cost
%             BestSoFar = BestCandidate;
%         end
    else
        ConsecFails = ConsecFails + 1;
        if ConsecFails > RestartCount
            % Restart
            ConsecFails = 0;
%             Population(1) = BestSoFar;
            T0 = T0*(-psi(1)); %Euler-Mascheroni constant 
            T = T0;
%             if (GenIndex > deriv)
%                 derivation = diff(MinCost(GenIndex-deriv:GenIndex));
%                 if abs(sum(derivation)) < eps
%                     TGama = TGama*0.5;
%                     T0 = TGama;
%                 end
%             end
%        elseif rand < exp((Population(1).cost - Candidate(1).cost) / T)
        elseif rand < exp(-1 / T)
            j = length(Population);
            for i = 1:round(cross * length(Candidate))
                if Population(j).cost > Candidate(i).cost
                    Population(j) = Candidate(i);
                    j = j - 1;
                else
                    break;
                end
            end
            Population = SortByCost(Population);
        end
    end
    DisplayThisTime = DisplayFlag && (mod(GenIndex, 100) == 0);
    TArray(GenIndex+1) = T;
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayThisTime);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
if DisplayFlag
    xlabel('Iteration')
    ylabel('Best Cost So Far')
    figure
    plot(0:GenIndex, TArray)
    xlabel('Iteration'), ylabel('Temperature')
end
return
