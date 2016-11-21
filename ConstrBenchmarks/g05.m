function [InitFunction, CostFunction, FeasibleFunction] = g05
% g05 Function
% Dimensions: 2
% Domain: Variable
% min x: 679.9453, 1026.067, 0.1188764, -0.3962336
% min f(x): 5126.4981 best known
% Constraints: 5
InitFunction = @g05Init;
CostFunction = @g05Cost;
FeasibleFunction = [];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g05Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue= 5126.4981;
OPTIONS.numVar = 4; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [0 0 -0.55 -0.55];
OPTIONS.MaxDomain = [1200 1200 0.55 0.55];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g05Cost(Population, ~)
popsize = length(Population);
epsilon = 0.0001; %inequality epsilon
d1 = 0.000002 / 3;
for popindex = 1 : popsize
    Population(popindex).cost = 3*Population(popindex).chrom(1) + ...
        0.000001*(Population(popindex).chrom(1))^3 + ...
        2*Population(popindex).chrom(2) + ...
        d1*(Population(popindex).chrom(2))^3;
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = -Population(popindex).chrom(4) + Population(popindex).chrom(3) - 0.55;
    %g2 <= 0
    Population(popindex).g(2) = -Population(popindex).chrom(3) + Population(popindex).chrom(4) - 0.55;
    %g3 <= 0, g4 <= 0
    %Common Method: Convert equality to 2 inequalities
    h3 = 1000*sin(-Population(popindex).chrom(3)-0.25) + ...
        1000*sin(-Population(popindex).chrom(4)-0.25)+894.8 - Population(popindex).chrom(1);
    Population(popindex).g(3) = h3 - epsilon;
    Population(popindex).g(4) = -h3 - epsilon;
    %g5 <= 0, g6 <= 0 
    h4 = 1000*sin(Population(popindex).chrom(3)-0.25)+...
        1000*sin(Population(popindex).chrom(3)-Population(popindex).chrom(4)-0.25)+894.8-...
        Population(popindex).chrom(2);
    Population(popindex).g(5) = h4 - epsilon;
    Population(popindex).g(6) = -h4 - epsilon;    
    %g7 <= 0, g8 <= 0
    h5 = 1000*sin(Population(popindex).chrom(4)-0.25)+...
        1000*sin(Population(popindex).chrom(4)-Population(popindex).chrom(3)-0.25)...
        +1294.8;
    Population(popindex).g(7) = h5 - epsilon;
    Population(popindex).g(8) = -h5 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return