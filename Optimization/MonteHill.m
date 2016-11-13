function [MeanMin, MeanMinNorm, BestMin, BestMinNorm] = MonteHill
% Monte Carlo simulations of hill climbing software
% OUTPUT MeanMin is the mean of the best solution found. It is an
%        nFunction x nBench array, where nFunction is the number of optimization
%        functions that are used, and nBench is the number of benchmarks that
%        are optimized.
% OUTPUT MeanMinNorm is MeanMin normalized to a minimum of 1 for each benchmark.
% OUTPUT BestMin is the best solution found by each optimization function
%        for each benchmark.
% OUTPUT BestMinNorm is BestMin normalized to a minimum of 1 for each benchmark.
nMonte = 50; % number of Monte Carlo runs
% Benchmark functions
 Bench = [     %     multimodal? separable?  regular?
 'Ackley     '; %     y           n           y
 'Fletcher   '; %     y           n           n
 'Griewank   '; %     y           n           y
 'Penalty1   '; %     y           n           y
 'Penalty2   '; %     y           n           y
 'Quartic    '; %     n           y           y
 'Rastrigin  '; %     y           y           y
 'Rosenbrock '; %     n           n           y
 'Schwefel12 '; %     y           y           n
 'Schwefel221'; %     n           n           y
 'Schwefel222'; %     y           n           n
 'Schwefel226'; %     n           n           n
 'Sphere     '; %     n           y           y
 'Step       ']; %    n           y           n
nBench = size(Bench, 1);
MeanMin = zeros(4, nBench);
BestMin = inf(4, nBench);
for j = 1 : nBench
    disp(['Benchmark function ', num2str(j), '/', num2str(nBench)]);
    for k = 1 : nMonte
        % Steepest ascent hill climbing
        [Cost] = eval(['SteepestHillClimbing(@', Bench(j,:), ', false);']);
        MeanMin(1,j) = ((k - 1) * MeanMin(1,j) + min(Cost)) / k;
        BestMin(1,j) = min(BestMin(1,j), min(Cost));
        % Next ascent ascent hill climbing
        [Cost] = eval(['NextHillClimbing(@', Bench(j,:), ', false);']);
        MeanMin(2,j) = ((k - 1) * MeanMin(2,j) + min(Cost)) / k;
        BestMin(2,j) = min(BestMin(2,j), min(Cost));
        % Random mutation hill climbing
        [Cost] = eval(['RandomHillClimbing(@', Bench(j,:), ', false);']);
        MeanMin(3,j) = ((k - 1) * MeanMin(3,j) + min(Cost)) / k;
        BestMin(3,j) = min(BestMin(3,j), min(Cost));
        % Adaptive hill climbing
        [Cost] = eval(['AdaptiveHillClimbing(@', Bench(j,:), ', false);']);
        MeanMin(4,j) = ((k - 1) * MeanMin(4,j) + min(Cost)) / k;
        BestMin(4,j) = min(BestMin(4,j), min(Cost));
    end
end
% Normalize the results
if min(MeanMin) == 0
    MeanMinNorm = [];
else
    MeanMinNorm = MeanMin * diag(1./min(MeanMin));
end
if min(BestMin) == 0
    BestMinNorm = [];
else
    BestMinNorm = BestMin * diag(1./min(BestMin));
end