function [InitFunction, CostFunction, FeasibleFunction] = g22
% g22 Function
% Dimensions: 22
% Domain: Variable
% min x: 236.430975504001054, 135.82847151732463, 204.818152544824585,
% 6446.54654059436416, 3007540.83940215595, 4074188.65771341929,
% 32918270.5028952882, 130.075408394314167, 170.817294970528621,
% 299.924591605478554, 399.258113423595205, 330.817294971142758,
% 184.51831230897065, 248.64670239647424, 127.658546694545862,
% 269.182627528746707, 160.000016724090955, 5.29788288102680571,
% 5.13529735903945728, 5.59531526444068827, 5.43444479314453499, 5.07517453535834395
% min f(x): 236.430975504001
% Constraints: 20
InitFunction = @g22Init;
CostFunction = @g22Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g22Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 236.430975504001;
OPTIONS.numVar = 22; % Problem dimension
%Boundaries for each variable
OPTIONS.MinDomain = [    0  0    0    0      0      0      0   100    100    100.01 100 100   0   0   0   0.01   0.01 -4.7  -4.7  -4.7  -4.7  -4.7];
OPTIONS.MaxDomain = [20000 10^6 10^6 10^6 4*10^7 4*10^7 4*10^7 299.99 399.99 300    400 600 500 500 500 300    400     6.25  6.25  6.25  6.25  6.25];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g22Cost(Population, ~)
popsize = length(Population);
epsilon = 0.0001; %inequality epsilon
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    Population(popindex).cost = x(1);
    % Compute the level of constraint violation
    Population(popindex).g(1) = -x(1) + power(x(2),0.6) + power(x(3),0.6) + power(x(4),0.6);
    %Common Method: Convert equality to 2 inequalities
    % h1 = 0
    h1 = x(5) - 100000*x(8) + 1*power(10,7);
    Population(popindex).g(2) = h1 - epsilon;
    Population(popindex).g(3) = -h1 - epsilon;
    % h2 = 0
    h2 = x(6) + 100000*x(8) - 100000*x(9);
    Population(popindex).g(4) = h2 - epsilon;
    Population(popindex).g(5) = -h2 - epsilon;    
    % h3 = 0
    h3 = x(7) + 100000*x(9) - 5*power(10,7);
    Population(popindex).g(6) = h3 - epsilon;
    Population(popindex).g(7) = -h3 - epsilon;
    % h4 = 0
    h4 = x(5) + 100000*x(10) - 3.3*power(10,7);
    Population(popindex).g(8) = h4 - epsilon;
    Population(popindex).g(9) = -h4 - epsilon;
    % h5 = 0
    h5 = x(6) + 100000*x(11) - 4.4*power(10,7);
    Population(popindex).g(10) = h5 - epsilon;
    Population(popindex).g(11) = -h5 - epsilon;    
    % h6 = 0
    h6 = x(7) + 100000*x(12) - 6.6*power(10,7);
    Population(popindex).g(12) = h6 - epsilon;
    Population(popindex).g(13) = -h6 - epsilon;
    % h7 = 0
    h7 = x(5) - 120*x(2)*x(13);
    Population(popindex).g(14) = h7 - epsilon;
    Population(popindex).g(15) = -h7 - epsilon;    
    % h8 = 0
    h8 = x(6) - 80*x(3)*x(14);
    Population(popindex).g(16) = h8 - epsilon;
    Population(popindex).g(17) = -h8 - epsilon;
    % h9 = 0
    h9 = x(7) - 40*x(4)*x(15);
    Population(popindex).g(18) = h9 - epsilon;
    Population(popindex).g(19) = -h9 - epsilon;
    % h10 = 0
    h10 = x(8) - x(11) + x(16);
    Population(popindex).g(20) = h10 - epsilon;
    Population(popindex).g(21) = -h10 - epsilon;    
    % h11 = 0
    h11 = x(9) - x(12) + x(17);
    Population(popindex).g(22) = h11 - epsilon;
    Population(popindex).g(23) = -h11 - epsilon;
    % h12 = 0
    h12 = -x(18) + log(x(10) - 100);
    Population(popindex).g(24) = h12 - epsilon;
    Population(popindex).g(25) = -h12 - epsilon;    
    % h13 = 0
    h13 = -x(19) + log(-x(8) + 300);
    Population(popindex).g(26) = h13 - epsilon;
    Population(popindex).g(27) = -h13 - epsilon;
    % h14 = 0
    h14 = -x(20) + log(x(16));
    Population(popindex).g(28) = h14 - epsilon;
    Population(popindex).g(29) = -h14 - epsilon;
    % h15 = 0
    h15 = -x(21) + log(-x(9) + 400);
    Population(popindex).g(30) = h15 - epsilon;
    Population(popindex).g(31) = -h15 - epsilon;    
    % h16 = 0
    h16 = -x(22) + log(x(17));
    Population(popindex).g(32) = h16 - epsilon;
    Population(popindex).g(33) = -h16 - epsilon;
    % h17 = 0
    h17 = -x(8) - x(10) + x(13)*x(18) - x(13)*x(19) + 400;
    Population(popindex).g(34) = h17 - epsilon;
    Population(popindex).g(35) = -h17 - epsilon;    
    % h18 = 0
    h18 = x(8) - x(9) - x(11) + x(14)*x(20) - x(14)*x(21) + 400;
    Population(popindex).g(36) = h18 - epsilon;
    Population(popindex).g(37) = -h18 - epsilon;
    % h19 = 0
    h19 = x(9) - x(12) - 4.60517*x(15) + x(15)*x(22) + 100;
    Population(popindex).g(38) = h19 - epsilon;
    Population(popindex).g(39) = -h19 - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
end
return