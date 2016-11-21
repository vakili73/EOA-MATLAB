function [InitFunction, CostFunction, FeasibleFunction] = g16
% g16 Function 
% Dimensions: 5
% Domain: Variable
% min x: 705:174537070090537,68.5999999999999943, 102.899999999999991,
% 282.324931593660324, 37.5841164258054832
% min f(x): -1.90515525853479
% Number of constraints: 38
InitFunction = @g16Init;
CostFunction = @g16Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = g16Init(OPTIONS)
% Initialize population
OPTIONS.MinFunValue = -1.90515525853479;
OPTIONS.numVar = 5; % Problem dimension
% Boundaries for each variable
OPTIONS.MinDomain = [704.4148  68.6    0    193      25];
OPTIONS.MaxDomain = [906.3855 288.88 134.75 287.0966 84.1988];
Population(OPTIONS.popsize).chrom = zeros(1, OPTIONS.numVar);
for popindex = 1 : OPTIONS.popsize
	chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1, OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
OPTIONS.OrderDependent = true;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = g16Cost(Population, ~)
% Compute the cost of each member in Population
popsize = length(Population);
for popindex = 1 : popsize
    x = Population(popindex).chrom;
    y = [];
    c = [];
    y(1) = x(2) + x(3) + 41.6;
    c(1) = 0.024*x(4) - 4.62;
    y(2) = 12.5/c(1) + 12;
    c(2) = 0.0003535*(x(1)^2) + 0.5311*x(1) + 0.08705*y(2)*x(1);
    c(3) = 0.052*x(1) + 78 + 0.002377*y(2)*x(1);
    y(3) = c(2)/c(3);
    y(4) = 19*y(3);
    c(4) = 0.04782*(x(1)-y(3)) + (0.1956*((x(1)-y(3))^2))/x(2) + 0.6376*y(4) + 1.594*y(3);
    c(5) = 100*x(2);
    c(6) = x(1) - y(3) - y(4);
    c(7) = 0.950 - c(4)/c(5);
    y(5) = c(6)*c(7);
    y(6) = x(1) - y(5) - y(4) - y(3);
    c(8) = (y(5) + y(4))*0.995;
    y(7) = c(8)/y(1);
    y(8) = c(8)/3798;
    c(9) = y(7) - (0.0663*y(7))/y(8) - 0.3153;
    y(9) = (96.82/c(9)) + 0.321*y(1);
    y(10) = 1.29*y(5) + 1.258*y(4) + 2.29*y(3) +1.71*y(6);
    y(11) = 1.71*x(1) - 0.452*y(4) + 0.580*y(3);
    c(10) = 12.3/752.3;
    c(11) = (1.75*y(2))*(0.995*x(1));
    c(12) = 0.995*y(10) + 1998;
    y(12) = c(10)*x(1) + (c(11)/c(12));
    y(13) = c(12) - 1.75*y(2);
    y(14) = 3623 + 64.4*x(2) + 58.4*x(3) + 146312/(y(9)+x(5));
    c(13) = 0.995*y(10) + 60.8*x(2) + 48*x(4) - 0.1121*y(14) - 5095;
    y(15) = y(13)/c(13);
    y(16) = 148000 - 331000*y(15) + 40*y(13) - 61*y(15)*y(13);
    c(14) = 2324*y(10) - 28740000*y(2);
    y(17) = 14130000 - 1328*y(10) - 531*y(11) + c(14)/c(12);
    c(15) = y(13)/y(15) - y(13)/0.52;
    c(16) = 1.104 - 0.72*y(15);
    c(17) = y(9) + x(5);
    Population(popindex).cost = 0.000117*y(14) + 0.1365 + 0.000001502*y(16) + 0.0321*y(12)....
        + 0.004324*y(5) + 0.0001*(c(15)/c(16)) + (37.48*y(2))/c(12) - 0.0000005843*y(17);
    % Compute the level of constraint violation
    %g1 <= 0
    Population(popindex).g(1) = (0.28/0.72)*y(5) - y(4);
    %g2 <= 0
    Population(popindex).g(2) = x(3) - 1.5*x(2);
    %g3 <= 0
    Population(popindex).g(3) = 3496*(y(2)/c(12)) - 21;
    %g4 <= 0
    Population(popindex).g(4) = 110.6 + y(1) - 62212/c(17);
    %g5 <= 0
    Population(popindex).g(5) = 213.1 - y(1);
    %g6 <= 0
    Population(popindex).g(6) = y(1) - 405.23;
    %g7 <= 0
    Population(popindex).g(7) = 17.505 - y(2);
    %g8 <= 0
    Population(popindex).g(8) = y(2) - 1053.6667;
    %g9 <= 0
    Population(popindex).g(9) = 11.275 - y(3);
    %g10 <= 0
    Population(popindex).g(10) = y(3) - 35.03;
    %g11 <= 0
    Population(popindex).g(11) = 214.228 - y(4);
    %g12 <= 0
    Population(popindex).g(12) = y(4) - 665.585;
    %g13 <= 0
    Population(popindex).g(13) = 7.458 - y(5);
    %g14 <= 0
    Population(popindex).g(14) = y(5) - 584.463;
    %g15 <= 0
    Population(popindex).g(15) = 0.961 - y(6);
    %g16 <= 0
    Population(popindex).g(16) = y(6) - 265.916;
    %g17 <= 0
    Population(popindex).g(17) = 1.612 - y(7);
    %g18 <= 0
    Population(popindex).g(18) = y(7) - 7.046;
    %g19 <= 0
    Population(popindex).g(19) = 0.146 - y(8);
    %g20 <= 0
    Population(popindex).g(20) = y(8) - 0.222;
    %g21 <= 0
    Population(popindex).g(21) = 107.99 - y(9);
    %g22 <= 0
    Population(popindex).g(22) = y(9) - 273.366;
    %g23 <= 0
    Population(popindex).g(23) = 922.693 - y(10);
    %g24 <= 0
    Population(popindex).g(24) = y(10) - 1286.105;
    %g25 <= 0
    Population(popindex).g(25) = 926.832 - y(11);
    %g26 <= 0
    Population(popindex).g(26) = y(11) - 1444.046;
    %g27 <= 0
    Population(popindex).g(27) = 18.766 - y(12);
    %g28 <= 0
    Population(popindex).g(28) = y(12) - 537.141;
    %g29 <= 0
    Population(popindex).g(29) = 1072.163 - y(13);
    %g30 <= 0
    Population(popindex).g(30) = y(13) - 3247.039;
    %g31 <= 0
    Population(popindex).g(31) = 8961.448 - y(14);
    %g32 <= 0
    Population(popindex).g(32) = y(14) - 26844.086;
    %g33 <= 0
    Population(popindex).g(33) = 0.063 - y(15);
    %g34 <= 0
    Population(popindex).g(34) = y(15) - 0.386;
    %g35 <= 0
    Population(popindex).g(35) = 71084.33 - y(16);
    %g36 <= 0
    Population(popindex).g(36) = -140000 + y(16);
    %g37 <= 0
    Population(popindex).g(37) = 2802713 - y(17);
    %g38 <= 0
    Population(popindex).g(38) = y(17) - 12146108;
    Population(popindex).g(Population(popindex).g <= 0) = 0;
    %Sum of constraints: Bigger the sum, bigger the constraint violation
    Population(popindex).G = sum(Population(popindex).g);
    %Geometric Mean of constrains: Bigger the value, bigger the constrain violation
end
return