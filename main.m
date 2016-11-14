clc;
clear;
close all;

%% Initializing
% use this function for generate random number from clock puls:
% fix(sum(100*clock));
Init__;

%% Start the program

% chaotic map functions
Chaotic = {         % Method of generate random function
    'rand       '   ,'rand()'               ;
    'randn      '   ,'randn()'              ;
    'Bernoulli  '   ,'Bernoulli()'          ;
    'Chebyshev  '   ,'Chebyshev()'          ;
    'Circle     '   ,'Circle()'             ;
    'Gauss      '   ,'Gauss()'              ;
    'ICMIC      '   ,'ICMIC()'              ;
    'Logistic   '   ,'Logistic()'           ;
    'Piecewise  '   ,'Piecewise()'          ;
    'Singer     '   ,'Singer()'             ;
    'Sinus      '   ,'Sinus()'              ;
    'Sinusoidal '   ,'Sinusoidal()'         ;
    'Tent       '   ,'Tent()'               };
nChaotic = size(Chaotic, 1);
Results = zeros(nChaotic, 4);
parfor i = 1 : nChaotic
    Results(i ,:) = MonteSA(Chaotic{i, 2});
end

function [MeanMin, MeanMinNorm, BestMin, BestMinNorm] = MonteSA(Random)
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
MeanMin = zeros(1, nBench);
BestMin = inf(1, nBench);
for j = 1 : nBench
    disp(['Benchmark function ', num2str(j), '/', num2str(nBench), ' Chaotic Map Function : ', Random]);
    for k = 1 : nMonte
        % Simulated Annealing
        [Cost] = eval(['SA(@', Bench(j,:), ', false, ''', Random, ''');']);
        MeanMin(1,j) = ((k - 1) * MeanMin(1,j) + min(Cost)) / k;
        BestMin(1,j) = min(BestMin(1,j), min(Cost));
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

end