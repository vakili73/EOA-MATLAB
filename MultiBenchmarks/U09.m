function [InitFunction, CostFunction, FeasibleFunction] = U09
% Adapated from cec09.m, which was downloaded from www.ntu.edu.sg/home/EPNSugan/
InitFunction = @Init;
CostFunction = @Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
% Initialize the population
if ~isfield(OPTIONS, 'MinDomain')
    OPTIONS.MinDomain = [0, -1 * ones(1, OPTIONS.numVar-1)];
end
if ~isfield(OPTIONS, 'MaxDomain')
    OPTIONS.MaxDomain = ones(1, OPTIONS.numVar);
end
OPTIONS.OrderDependent = true; % the order of the independent variables does matter
OPTIONS = ComputeRandomShift(OPTIONS);
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
num = length(Population);
dim =  OPTIONS.numVar;
epsi = 0.1;
for pop = 1 : num
    x = Population(pop).chrom;
    % f1
    temp1 = 0;
    for j = 4 : 3 : dim
        temp1 = temp1 + (x(j) - ((2*x(2)*sin((6*pi*x(1)) + (j*pi/dim)))))^2;
    end
    J1 = 4 : 3 : dim;
    Population(pop).cost(1) = (0.5*(max(0,((1+epsi)*(1-(4*((2*x(1)-1)^2))))) + (2*x(1)))*x(2)) + (2 * temp1 / length(J1));
    % f2
    temp2 = 0;
    for j = 5 : 3 : dim
        temp2 = temp2 + (x(j) - ((2*x(2)*sin((6*pi*x(1)) + (j*pi/dim)))))^2;
    end
    J2 = 5 : 2 : dim;
    Population(pop).cost(2) = (0.5*(max(0,((1+epsi)*(1-(4*((2*x(1)-1)^2))))) - (2*x(1)) + 2)*x(2)) + (2 * temp2 / length(J2));
    % f3
    temp3 = 0;
    for j = 3 : 3 : dim
        temp3 = temp3 + (x(j) - ((2*x(2)*sin((6*pi*x(1)) + (j*pi/dim)))))^2;
    end
    J3 = 3 : 3 : dim;
    Population(pop).cost(3) = 1 - x(2) + (2 * temp3 / length(J3));
end
return