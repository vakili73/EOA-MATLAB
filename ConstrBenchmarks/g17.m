function [InitFunction, CostFunction, FeasibleFunction] = g17
% g17 Function
% Dimensions: 6
% Domain: Variable
% min x: 201.784467214523659, 99.9999999999999005, 383.071034852773266, 420, -10.9076584514292652, 0.0731482312084287128
% min f(x): 8853.53967480648
% Constraints: 4
InitFunction = @g17Init;
CostFunction = @g17Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g17Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 8853.53967480648;
OPTIONS.numVar = 6; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [  0    0 340 340 -1000 0];
OPTIONS.MaxDomain = [400 1000 420 420 +1000 0.5236];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g17Cost(Population, ~)
popsize = length(Population);
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    if 0 <= x(1) < 300
        f1 = 30*x(1);
    end
    if 300 <= x(1) < 400
        f1 = 31*x(1);
    end
    if 0 <= x(2) < 100
        f2 = 28*x(2);
    end
    if 100 <= x(1) < 200
        f2 = 30*x(1);
    end
    if 200 <= x(1) < 1000
        f2 = 30*x(1);
    end
    Population(popindex).cost = f1 + f2;
    % Compute the level of constraint violation
    %Common Method: Convert equality to 2 inequalities
    % h1 = 0
    h1 = -x(1) + 300 - ((x(3)*x(4))/131.078)*cos(1.48477 - x(6)) + ((0.90798*power(x(3),2))/131.078)*cos(1.47588);
    Population(popindex).g(1) = h1 - epsilon;
    Population(popindex).g(2) = -h1 - epsilon;
    % h2 = 0
    h2 = -x(2) - ((x(3)*x(4))/131.078)*cos(1.48477 + x(6)) + ((0.90798*power(x(4),2))/131.078)*cos(1.47588);
    Population(popindex).g(3) = h2 - epsilon;
    Population(popindex).g(4) = -h2 - epsilon;    
    % h3 = 0
    h3 = -x(5) - ((x(3)*x(4))/131.078)*sin(1.48477 + x(6)) + ((0.90798*power(x(3),2))/131.078)*sin(1.47588);
    Population(popindex).g(5) = h3 - epsilon;
    Population(popindex).g(6) = -h3 - epsilon;
    % h4 = 0
    h4 = 200 - ((x(3)*x(4))/131.078)*sin(1.48477 - x(6)) + ((0.90798*power(x(3),2))/131.078)*sin(1.47588);
    Population(popindex).g(7) = h4 - epsilon;
    Population(popindex).g(8) = -h4 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return