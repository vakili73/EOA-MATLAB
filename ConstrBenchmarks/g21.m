function [InitFunction, CostFunction, FeasibleFunction] = g21
% g21 Function
% Dimensions: 7
% Domain: Variable
% min x: 193.724510070034967,5.56944131553368433e-27,17.3191887294084914,100.04789
%       7801386839,6.68445185362377892,5.99168428444264833,6.21451648886070451
% min f(x): 193.724510070035
% Constraints: 6
InitFunction = @g21Init;
CostFunction = @g21Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g21Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 193.724510070035;
OPTIONS.numVar = 7; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [   0  0  0 100 6.3 5.9 4.5];
OPTIONS.MaxDomain = [1000 40 40 300 6.7 6.4 6.25];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g21Cost(Population, ~)
popsize = length(Population);
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = x(1);
    % Compute the level of constraint violation
    Population(popindex).g(1) = -x(1) + 35*power(x(2),0.6) + 35*power(x(3),0.6);
    %Common Method: Convert equality to 2 inequalities
    % h1 = 0
    h1 = -300*x(3) + 7500*x(5) - 7500*x(6) - 25*x(4)*x(5) + 25*x(4)*x(6) + x(3)*x(4);
    Population(popindex).g(2) = h1 - epsilon;
    Population(popindex).g(3) = -h1 - epsilon;
    % h2 = 0
    h2 = 100*x(2) + 155.365*x(4) + 2500*x(7) - x(2)*x(4) - 25*x(4)*x(7) - 15536.5;
    Population(popindex).g(4) = h2 - epsilon;
    Population(popindex).g(5) = -h2 - epsilon;    
    % h3 = 0
    h3 = -x(5) + log(-x(4) + 900);
    Population(popindex).g(6) = h3 - epsilon;
    Population(popindex).g(7) = -h3 - epsilon;
    % h4 = 0
    h4 = -x(6) + log(x(4) + 300);
    Population(popindex).g(8) = h4 - epsilon;
    Population(popindex).g(9) = -h4 - epsilon;
    % h5 = 0
    h5 = -x(7) + log(-2*x(4) + 700);
    Population(popindex).g(10) = h5 - epsilon;
    Population(popindex).g(11) = -h5 - epsilon;    
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return