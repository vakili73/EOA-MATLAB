function Init__
% add all requered path to matlab
if isempty(strfind(path, ['Benchmarks', pathsep]))
    addpath('Benchmarks')
end
if isempty(strfind(path, ['BenchPlot', pathsep]))
    addpath('BenchPlot')
end
if isempty(strfind(path, ['ChaoticMap', pathsep]))
    addpath('ChaoticMap')
end
if isempty(strfind(path, ['ChaoticPlot', pathsep]))
    addpath('ChaoticPlot')
end
if isempty(strfind(path, ['CommonCode', pathsep]))
    addpath('CommonCode')
end
if isempty(strfind(path, ['GA', pathsep]))
    addpath('GA')
end
if isempty(strfind(path, ['GAMath', pathsep]))
    addpath('GAMath')
end
if isempty(strfind(path, ['PSO', pathsep]))
    addpath('PSO')
end
if isempty(strfind(path, ['Optimization', pathsep]))
    addpath('Optimization')
end
if isempty(strfind(path, ['SimulatedAnnealing', pathsep]))
    addpath('SimulatedAnnealing')
end
if isempty(strfind(path, ['ConstrBenchmarks', pathsep]))
    addpath('ConstrBenchmarks')
end
if isempty(strfind(path, ['MultiBenchmarks', pathsep]))
    addpath('MultiBenchmarks')
end

end

