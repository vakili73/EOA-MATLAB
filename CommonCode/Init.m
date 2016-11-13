function [OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed)
% Initialize population-based optimization software.
% WARNING: Some of the optimization routines may not work if population size is odd.
if ~exist('OPTIONS', 'var') || isempty(OPTIONS) || ~isfield(OPTIONS, 'popsize')
    OPTIONS.popsize = 100; % total population size
end
if ~isfield(OPTIONS, 'Maxgen')
    OPTIONS.Maxgen = 100; % generation count limit
end
if ~isfield(OPTIONS, 'pmutate')
    OPTIONS.pmutate = 0.01; % mutation probability
end
if ~isfield(OPTIONS, 'clearDups')
    OPTIONS.clearDups = true; % whether or not to replace duplicate individuals with random individuals
end
if ~isfield(OPTIONS, 'Keep')
    OPTIONS.Keep = 2; % number of elites
end
if ~isfield(OPTIONS, 'numVar')
    OPTIONS.numVar = 10; % number of dimensions in cost function
end
if ~isfield(OPTIONS, 'ShiftFlag')
    OPTIONS.ShiftFlag = 0; % whether or not to randomly shift the solution of the optimization problem
end
if ~isfield(OPTIONS, 'RotationFlag')
    OPTIONS.RotationFlag = 0; % whether or not to randomly rotate the solution of the optimization problem
end
if ~isfield(OPTIONS, 'ConstrMethod') || isempty(OPTIONS.ConstrMethod)
%    OPTIONS.ConstrMethod = 1; % Deb's niched-penalty approach
%    OPTIONS.ConstrMethod = 2; % Eclectic GA
%    OPTIONS.ConstrMethod = 3; % Dynamic constraints combined with superiority of feasible points
%    OPTIONS.ConstrMethod = 4; % Dynamic constraints
%    OPTIONS.ConstrMethod = 5; % Exponential dynamic penalty combined with superiority of feasible points
%    OPTIONS.ConstrMethod = 6; % Exponential dynamic penalty
%    OPTIONS.ConstrMethod = 7; % Adaptive penalty weights
%    OPTIONS.ConstrMethod = 8; % Stochastic ranking
    OPTIONS.ConstrMethod = 9; % Epsilon-level comparisons
end
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = fix(sum(100*clock));
end
SetPlotOptions
% Initialize the random number generator. The initialization method depends on the Matlab version.
VersionStr = ver;
for i = 1 : length(VersionStr)
    if isequal(VersionStr(i).Name, 'MATLAB')
        StringVersion = VersionStr(i).Version;
        DecimalPosition = strfind(StringVersion, '.');
        IntegerVersion = str2double(StringVersion(1:DecimalPosition-1));
        FractionVersion = str2double(StringVersion(DecimalPosition+1:end));
        if (IntegerVersion < 7) || ((IntegerVersion == 7) && (FractionVersion < 7))
            rand('state', RandSeed); %#ok<RAND>
            randn('state', RandSeed); %#ok<RAND>
        else
            rng(RandSeed);
        end
        if DisplayFlag
            display(['RandSeed = ', num2str(RandSeed)]);
        end
        break
    end
end
% Get the addresses of the initialization, cost, and feasibility functions.
[OPTIONS.InitFunction, OPTIONS.CostFunction] = ProblemFunction();
% Initialize the population.
[Population, OPTIONS] = OPTIONS.InitFunction(OPTIONS);
% Compute cost of each individual  
Population = OPTIONS.CostFunction(Population, OPTIONS);
% Sort the population from most fit to least fit
Population = PopSort(Population, OPTIONS.ConstrMethod, 1);
MinCost = zeros(OPTIONS.Maxgen+1, 1);
AvgCost =  zeros(OPTIONS.Maxgen+1, 1);
% Keep track of the minimum and average number of constraint violations,
% but only for single-objective optimization problems
if length(Population(1).cost) > 1
    MinConstrViol = [];
    AvgConstrViol = [];
else    
    MinConstrViol = zeros(OPTIONS.Maxgen+1, 1);
    AvgConstrViol = zeros(OPTIONS.Maxgen+1, 1);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, 0, DisplayFlag);
end
return