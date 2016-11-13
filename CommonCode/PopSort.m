function [Population] = PopSort(Population, ConstrMethod, GenIndex)
% Sort the population members from best to worst
persistent epsilon0
if isfield(Population, 'G')
    % Constrained problem
    if ConstrMethod == 1
        % Deb's penalty approach (niched penalty approach)
        % Sort the feasible individuals by cost
        Feasibles = Population([Population.G] == 0);
        Feasibles = SortByCost(Feasibles);
        % Sort the infeasible individuals by the sum of their constraint violations
        Infeasibles = Population([Population.G] ~= 0);
        numInfeas = length(Infeasibles);
        numConstr = length(Population(1).g);
        if numInfeas > 0
            [G, indices] = sort([Infeasibles.G], 'ascend');
            Chroms = zeros(numInfeas, length(Infeasibles(1).chrom));
            Cost = zeros(numInfeas, 1);
            g = zeros(numInfeas, numConstr);
            for i = 1 : numInfeas
                Chroms(i, :) = Infeasibles(indices(i)).chrom;
                Cost(i) = Infeasibles(indices(i)).cost;
                g(i, :) = Infeasibles(indices(i)).g;
            end
            for i = 1 : numInfeas
                Infeasibles(i).G = G(i);
                Infeasibles(i).chrom = Chroms(i, :);
                Infeasibles(i).cost = Cost(i);
                Infeasibles(i).g = g(i, :);
            end
        end
        % The final sorted population is comprised of the sorted feasibles followed by the sorted infeasibles
        Population = [Feasibles, Infeasibles];
    elseif ConstrMethod == 2
        % Eclectic EA
        % Sort the feasible individuals by cost
        Feasibles = Population([Population.G] == 0);
        Feasibles = SortByCost(Feasibles);
        % Set the cost of each infeasible individual proportional to the number of constraint violations        
        Infeasibles = Population([Population.G] ~= 0);
        numInfeas = length(Infeasibles);
        if numInfeas > 0
            numConstr = length(Infeasibles(1).g);
            ViolMags = reshape([Infeasibles.g], numConstr, numInfeas)';
            CostOffset = 0;
            if ~isempty(Feasibles)
                CostOffset = max([Feasibles.cost]);
            end
            for i = 1 : numInfeas
                numSatisfy = sum(ViolMags(i, :) == 0);
                Infeasibles(i).cost = CostOffset + 1e6 * (numConstr - numSatisfy);
            end
            % Sort the infeasibles by their cost
            Infeasibles = SortByCost(Infeasibles);
        end
        % The final sorted population is comprised of the sorted feasibles followed by the sorted infeasibles
        Population = [Feasibles, Infeasibles];
    elseif ConstrMethod == 3
        % Dynamic penalty method combined with superiority of feasible points
        c_const = 10;
        alpha = 2;
        % Sort the feasible individuals by cost
        Feasibles = Population([Population.G] == 0);
        Feasibles = SortByCost(Feasibles);
        % Set the cost of each infeasible individual to a penalized cost function
        Infeasibles = Population([Population.G] ~= 0);
        Infeasibles = SortByDynamicPenalizedCost(Infeasibles, c_const, GenIndex, alpha);
        % The final sorted population is comprised of the sorted feasibles followed by the sorted infeasibles
        Population = [Feasibles, Infeasibles];
    elseif ConstrMethod == 4
        % Dynamic penalty method
        c_const = 10;
        alpha = 2;
        Population = SortByDynamicPenalizedCost(Population, c_const, GenIndex, alpha);
    elseif ConstrMethod == 5
        % Exponential dynamic penalty combined with superiority of feasible points
        alpha = 10;
        % Sort the feasible individuals by cost
        Feasibles = Population([Population.G] == 0);
        Feasibles = SortByCost(Feasibles);
        % Set the cost of each infeasible individual to a penalized cost function
        Infeasibles = Population([Population.G] ~= 0);
        Infeasibles = SortByExpPenalizedCost(Infeasibles, GenIndex, alpha);
        % The final sorted population is comprised of the sorted feasibles followed by the sorted infeasibles
        Population = [Feasibles, Infeasibles];
    elseif ConstrMethod == 6
        % Exponential dynamic penalty method
        alpha = 10;
        Population = SortByExpPenalizedCost(Population, GenIndex, alpha);
    elseif ConstrMethod == 7
        % Adaptive penalty weights - parameters taken from Hadj-Alouane, 1993 technical report
        beta1 = 4;
        beta2 = 3;
        k_const = length(Population(1).chrom);
        Population = SortByAdaptivePenalizedCost(Population, GenIndex, beta1, beta2, k_const);
    elseif ConstrMethod == 8
        % Stochastic ranking
        % Bubble sort modified to use stochastic ranking
        Pf = 0.45;
        PopSize = length(Population);
        while true
            newn = 0;
            for i = 1 : PopSize-1
                Swap1 = (Population(i).G == 0) && (Population(i+1).G == 0) && (Population(i).cost > Population(i+1).cost);
                Swap2 = (rand < Pf) && (Population(i).cost > Population(i+1).cost);
                Swap3 = ~Swap1 && ~Swap2 && (Population(i).G > Population(i+1).G);
                if Swap1 || Swap2 || Swap3
                    temp = Population(i);
                    Population(i) = Population(i+1);
                    Population(i+1) = temp;
                    newn = i + 1;
                end
            end
            PopSize = newn;
            if PopSize <= 1, break, end
        end
    elseif ConstrMethod == 9
        % Epsilon-level comparisons
        PopSize = length(Population);
        if GenIndex == 1
            ViolationMags = sort([Population.G], 'ascend');
            theta = round(0.2 * PopSize);
            epsilon0 = ViolationMags(theta);
        end
        cp = 100;
        Tc = 200;
        epsilon = epsilon0 * (1 - GenIndex / Tc)^cp;
        epsilon = max(epsilon, 1e-10);
        % Bubble sort modified to use epsilon-level comparisons on constraint violation magnitudes
        while true
            newn = 0;
            for i = 1 : PopSize-1
                Swap1 = (Population(i).cost > Population(i+1).cost) && (Population(i).G < epsilon) && (Population(i+1).G < epsilon);
                Swap2 = (Population(i).G > Population(i+1).G) && (Population(i).G > epsilon);
                if Swap1 || Swap2
                    temp = Population(i);
                    Population(i) = Population(i+1);
                    Population(i+1) = temp;
                    newn = i + 1;
                end
            end
            PopSize = newn;
            if PopSize <= 1, break, end
        end
    end
else
    % Unconstrained problem
    Population = SortByCost(Population);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SortByCost(Population)
% Sort Population by its cost field using a bubble sort algorithm
PopSize = length(Population);
if PopSize == 0, return, end
while true
    newn = 0;
    for i = 1 : PopSize-1
        if Population(i).cost > Population(i+1).cost
            temp = Population(i);
            Population(i) = Population(i+1);
            Population(i+1) = temp;
            newn = i + 1;
        end
    end
    PopSize = newn;
    if PopSize <= 1, break, end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = NormalizeConstrViol(Population)
% Normalize the sum of the constraint violations in the Infeasibles population to lie between 0 and 1
PopSize = length(Population);
ViolMax = max([Population.G]);
if ViolMax == 0, return, end
for i = 1 : PopSize
    Population(i).G = Population(i).G / ViolMax;
end
% The following code normalizes each individual constraint violation
% numConstr = length(Infeasibles(1).g);
% ViolMags = reshape([Infeasibles.g], numConstr, numInfeas)';
% ViolMax = max(ViolMags, [], 1);
% ViolMax(ViolMax == 0) = 1; % If any constraints have a max violation magnitude of 0, set them equal to 1
% ViolMags = ViolMags * diag(1 ./ ViolMax);
% for i = 1 : numInfeas
%     Infeasibles(i).g = ViolMags(i, :);
%     Infeasibles(i).G = sum(Infeasibles(i).g);
% end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = NormalizeCost(Population)
% Normalize the cost in the population to be between 0 and 1.
% Put the normalized cost in the costTemp field.
PopSize = length(Population);
CostMin = min([Population.cost]);
CostMax = max([Population.cost]);
if CostMax == CostMin
    for i = 1 : PopSize
        Population(i).costTemp = Population(i).cost;
    end
else
    for i = 1 : PopSize
        Population(i).costTemp = (Population(i).cost - CostMin) / (CostMax - CostMin);
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SortByDynamicPenalizedCost(Population, c_const, GenIndex, alpha)
PopSize = length(Population);
if PopSize > 0
    Population = NormalizeConstrViol(Population);
    Population = NormalizeCost(Population);
    for i = 1 : PopSize
        UnpenalizedCost = Population(i).cost; % save the unpenalized cost
        Population(i).cost = Population(i).costTemp + (c_const * GenIndex)^alpha * Population(i).G;
        Population(i).costTemp = UnpenalizedCost; % store the unpenalized cost in the costTemp field
    end
    % Sort the individuals by their penalized cost
    Population = SortByCost(Population);
    % Restore the unpenalized cost to the cost field, restore the unnormalized constraint violation sum
    for i = 1 : PopSize
        Population(i).cost = Population(i).costTemp;
        Population(i).G = sum(Population(i).g);
    end
    Population = rmfield(Population, 'costTemp');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SortByExpPenalizedCost(Population, GenIndex, alpha)
PopSize = length(Population);
if PopSize > 0
    Population = NormalizeConstrViol(Population);
    Population = NormalizeCost(Population);
    for i = 1 : PopSize
        UnpenalizedCost = Population(i).cost; % save the unpenalized cost
        Population(i).cost = Population(i).costTemp * exp(alpha * Population(i).G * sqrt(GenIndex));
        Population(i).costTemp = UnpenalizedCost; % store the unpenalized cost in the costTemp field
    end
    % Sort the individuals by their penalized cost
    Population = SortByCost(Population);
    % Restore the unpenalized cost to the cost field
    for i = 1 : PopSize
        Population(i).cost = Population(i).costTemp;
        Population(i).G = sum(Population(i).g);
    end
    Population = rmfield(Population, 'costTemp');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = SortByAdaptivePenalizedCost(Population, GenIndex, beta1, beta2, k_const)
persistent BestFeasCount NoFeasCount rweight
if GenIndex == 1
    BestFeasCount = 0; % number of consecutive generations that the lowest-cost individual was feasible
    NoFeasCount = 0; % number of consecutive generations that the lowest-cost individual was infeasible
    rweight = 1;
end
PopSize = length(Population);
if PopSize > 0
    [~, ndx] = min([Population.cost]);
    if Population(ndx).G == 0
        BestFeasCount = BestFeasCount + 1;
        NoFeasCount = 0;
        if BestFeasCount >= k_const
            rweight = rweight / beta1;
        end
    else
        BestFeasCount = 0;
        NoFeasCount = NoFeasCount + 1;
        if NoFeasCount >= k_const
            rweight = rweight * beta2;
        end
    end
    Population = NormalizeConstrViol(Population);
    Population = NormalizeCost(Population);
    for i = 1 : PopSize
        UnpenalizedCost = Population(i).cost; % save the unpenalized cost
        Population(i).cost = Population(i).costTemp + rweight * Population(i).G;
        Population(i).costTemp = UnpenalizedCost; % store the unpenalized cost in the costTemp field
    end
    % Sort the individuals by their penalized cost
    Population = SortByCost(Population);
    % Restore the unpenalized cost to the cost field
    for i = 1 : PopSize
        Population(i).cost = Population(i).costTemp;
        Population(i).G = sum(Population(i).g);
    end
    Population = rmfield(Population, 'costTemp');
end
return