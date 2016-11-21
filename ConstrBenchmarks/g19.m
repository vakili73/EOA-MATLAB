function [InitFunction, CostFunction, FeasibleFunction] = g19
% g19 Function 
% Dimensions: 15
% Domain: Variable
% min x: 1.66991341326291344e-17, 3.95378229282456509e-16, 3.94599045143233784, 1.060
%       36597479721211e-16, 3.2831773458454161, 9.99999999999999822, 1.1282941467160
%       5333e-17, 1.2026194599794709e-17, 2.50706276000769697e-15, 2.246241229879706
%       77e-15, 0.370764847417013987, 0.278456024942955571, 0.523838487672241171, 0.388620152510322781, 0.298156764974678579
% min f(x): 32.6555929502463
% Number of constraints: 5
InitFunction = @g19Init;
CostFunction = @g19Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g19Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 32.6555929502463;
OPTIONS.numVar = 15; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];
OPTIONS.MaxDomain = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g19Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
temp = 0;
temp2 = 0;
temp3 = 0;
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    a = [ -16   2   0   1    0;
            0  -2   0 0.4    2;
         -3.5   0   2   0    0;
            0  -2   0  -4   -1;
            0  -9  -2   1 -2.8;
            2   0  -4   0    0;
           -1  -1  -1  -1   -1;
           -1  -2  -3  -2   -1;
            1   2   3   4    5;
            1   1   1   1    1 ];
    b = [-40 -2 -0.25 -4 -4 -1 -40 -60 5 1];
    c = [  30 -20 -10  32  -10;
          -20  39  -6 -31   32;
          -10  -6  10  -6  -10;
           32 -31  -6  39  -20;
          -10  32 -10 -20   30 ];
    d = [   4   8  10   6    2 ];
    e = [ -15 -27 -36 -18  -12 ];
    for j = 1 : 5
        for i = 1 : 5
            temp = temp + c(i,j)*x(10+i)*x(10+j);
        end
    end
    for j = 1 : 5
        temp2 = temp2 + d(j)*(x(10+j)^3);
    end
    for i = 1 : 10
        temp3 = temp3 + b(i)*x(i);
    end
    Population(popindex).cost = temp + (2*temp2) - temp3;
    % Compute the level of constraint violation
    % g1 : g5 <= 0
    for j = 1 : 5
       sum1 = 0;
       for i = 1 : 5
         sum1 = sum1 + c(i,j) * x(10 + i);
       end
       sum2 = 0;
       for i = 1 : 10
         sum2 = sum2 + a(i,j) * x(i);
       end
       Population(popindex).g(j) = -(2.0 * sum1) - (3.0 * d(j) * (x(10 + j)^2)) - e(j) + sum2;
    end
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return