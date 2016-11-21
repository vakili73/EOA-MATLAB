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
    'rand       '   ,'rand()-rand()'        ;
    'randn      '   ,'randn()'              ;
    'Bernoulli  '   ,'Bernoulli()-Bernoulli()'          ;
    'Chebyshev  '   ,'Chebyshev()'          ;
    'Circle     '   ,'Circle()'             ;
%     'Gauss      '   ,'Gauss()'              ;
    'GaussMap   '   ,'GaussMap()-GaussMap()'           ;
%     'ICMIC      '   ,'ICMIC()'              ;
%     'Logistic   '   ,'Logistic()'           ;
%     'Piecewise  '   ,'Piecewise()'          ;
%     'Singer     '   ,'Singer()'             ;
%     'Sinus      '   ,'Sinus()'              ;
%     'Sinusoidal '   ,'Sinusoidal()'         ;
%     'Tent       '   ,'Tent()'               ;
    'CMFoCSA    '   ,'CMFoCSA()'            };
nChaotic = size(Chaotic, 1);
Results = cell(nChaotic, 1);
parfor i = 1 : nChaotic
    Results{i} = F(Chaotic{i, 2});
end

function Results = F(chaotic)
[MeanMin, MeanMinNorm, BestMin, BestMinNorm, STDs] = MonteSA(chaotic);
Results = cell(1, 5);
Results{1, 1} = MeanMin;
Results{1, 2} = MeanMinNorm;
Results{1, 3} = BestMin;
Results{1, 4} = BestMinNorm;
Results{1, 5} = STDs;
end

function [MeanMin, MeanMinNorm, BestMin, BestMinNorm, STDs] = MonteSA(Random)
% Monte Carlo simulations of Simulated Anealing software
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
    %     'Penalty1   '; %     y           n           y
    %     'Penalty2   '; %     y           n           y
    'Quartic    '; %     n           y           y
    'Rastrigin  '; %     y           y           y
    'Rosenbrock '; %     n           n           y
    %     'Schwefel12 '; %     y           y           n
    %     'Schwefel221'; %     n           n           y
    %     'Schwefel222'; %     y           n           n
    'Schwefel226'; %     n           n           n
    'Sphere     '; %     n           y           y
    'Step       ']; %    n           y           n
nBench = size(Bench, 1);
MeanMin = zeros(1, nBench);
BestMin = inf(1, nBench);
STDs = zeros(1, nBench);
for j = 1 : nBench
    MinCosts = zeros(1, nMonte);
    for k = 1 : nMonte
        % Simulated Annealing
        [Cost] = eval(['SA(@', Bench(j,:), ', false, ''', Random, ''');']);
        MeanMin(1,j) = ((k - 1) * MeanMin(1,j) + min(Cost)) / k;
        BestMin(1,j) = min(BestMin(1,j), min(Cost));
        MinCosts(1,k) = min(Cost);
    end
    STDs(1,j) = std(MinCosts);
    disp(['BencFunc ', num2str(j), '/', num2str(nBench), ': ',Bench(j,:), ' MeanMin: ', num2str(MeanMin(1,j)), ' STD: ', num2str(STDs(1,j)),' CMFunc : ', Random]);
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