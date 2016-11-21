function [Pops] = EnumPops(n, N, LoopLimit, Pops)

% Recursively enumerate all possible populations.
% Since this routine is recursive, if you pass too high of a value in n and/or N, you will run out of memory.
% INPUTS: n = cardinality of search space, i.e., number of unique individuals in search space
%         N = population size
%         LoopLimit = ... well, loop limit, of course - don't pass anything in this variable on the first call
%         Pops = population distributions which have been calculated so far
% OUTPUTS: Pops = two-dimensional array of all possible population distributions
%                 Each row corresponds to one population
%                 The number in each column corresonds to the quantity of that particular individual

global NumCalls index iArray
if ~exist('LoopLimit', 'var')
    % first call
    NumPossible = nchoosek(n+N-1, N); % total number of possible populations
    Pops = zeros(NumPossible, n); % Each row of Pops will contain a population distribution
    index = 0; % Pops row index that has most recently been written
    NumCalls = 0; % Number of recursive calls to this routine
    iArray = zeros(1, n); % Population distribution which will be written to Pops (if distribution is valid)
    Pops = EnumPops(n, N, N, Pops);
else
    NumCalls = NumCalls + 1;
    for i = 0 : LoopLimit
        iArray(NumCalls) = i;
        if (NumCalls == n)
            if (sum(iArray) == N)
                index = index + 1;
                Pops(index, :) = iArray;
            end
        else
            Pops = EnumPops(n, N, LoopLimit-i, Pops);
            NumCalls = NumCalls - 1;
        end
    end
end
return