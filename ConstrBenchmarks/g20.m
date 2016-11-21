function [InitFunction, CostFunction, FeasibleFunction] = g20
% g19 Function 
% Dimension: 24
% Domain: Variable
% min x: 1.28582343498528086e-18,4.83460302526130664e-34,0,0,
% 6.30459929660781851e-18,7.57192526201145068e-34,
% 5.03350698372840437e-34,9.28268079616618064e-34,0,1.76723384525547359e-17,3.55686101822965701e-34,2.99413850083471346e-34,0.158143376337580827,2.2
% 9601774161699833e-19,1.06106938611042947e-18,1.31968344319506391e-18,0.53
% 0902525044209539,0,2.89148310257773535e-18,3.34892126180666159e-18,0,0.31
% 0999974151577319,5.41244666317833561e-05,4.84993165246959553e-16
% min f(x): 32.6555929502463
% Number of constraints: 20
InitFunction = @g20Init;
CostFunction = @g20Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g20Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = 32.6555929502463;
OPTIONS.numVar = 24; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0];
OPTIONS.MaxDomain = [10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g20Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
temp = 0;
epsilon = 0.0001;
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    a = [ 0.0693 0.0577 0.05 0.2 0.26 0.55 0.06 0.1 0.12 0.18 0.1 0.09 0.0693 0.0577 0.05 0.2 0.26 0.55 0.06 0.1 0.12 0.18 0.1 0.09];
    b = [ 44.094 58.12 58.12 137.4 120.9 170.9 62.501 84.94 133.425 82.507 46.07 60.097 
        44.094 58.12 58.12 137.4 120.9 170.9 62.501 84.94 133.425 82.507 46.07 60.097];
    c = [ 123.7 31.7 45.7 14.7 84.7 27.7 49.7 7.1 2.1 17.7 0.85 0.64];
    d = [ 31.244 36.12 34.784 92.7 82.7 91.6 56.708 82.7 80.8 64.517 49.4 49.1];
    e = [ 0.1 0.3 0.4 0.3 0.6 0.3];
    for i = 1 : 24 
        temp = temp + a(i)*x(i);
    end    
    Population(popindex).cost = temp;
    % Compute the level of constraint violation
    % g1 : g5 <= 0
    for i = 1 : 6
       sum1 = 0;
       for j = 1 : 24
         sum1 = sum1 + x(j) + e(i);
       end
       if i < 4
           Population(popindex).g(j) = (x(i) + x(i+12))/sum1;
       else
           Population(popindex).g(j) = (x(i+3) + x(i+12))/sum1;
       end
    end
    k = 7;
    h = zeros(1, 12);
    for i = 1 : 12
        sum1 = 0;
        sum2 = 0;
        for j = 1 : 24
            if j < 13
                sum1 = sum1 + (x(j)/b(j));
            else
                sum2 = sum2 + (x(j)/b(j));
            end
        end
        h(i) = x(i+12)/(b(i+12)*sum2) - c(i)*x(i)/(40*b(i)*sum1);
        Population(popindex).g(k) = h(i) - epsilon;
        k = k + 1;
        Population(popindex).g(k) = -h(i) - epsilon;
        k = k + 1;
    end
    sum3 = 0;
    for i = 1 : 24
        sum3 = sum3 + x(i);
    end
    h(13) = sum3 - 1;
    Population(popindex).g(31) = h(13) - epsilon;
    Population(popindex).g(32) = -h(13) - epsilon;
    sum4 = 0;
    sum5 = 0;
    for j = 1 : 24
        if j < 13
            sum4 = sum4 + (x(j)/d(j));
        else
            sum5 = sum5 + (x(j)/b(j));
        end
    end
    h(14) = sum4 + k*sum5 - 1.671;
    Population(popindex).g(33) = h(14) - epsilon;
    Population(popindex).g(34) = -h(14) - epsilon;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return