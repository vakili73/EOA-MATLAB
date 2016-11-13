function [InitFunction, CostFunction, FeasibleFunction] = AckleyDisc
InitFunction = @Init;
CostFunction = @Cost;
FeasibleFunction = @Feasible;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
if ~isfield(OPTIONS, 'Dim')
    OPTIONS.Dim = 2;
end
if ~isfield(OPTIONS, 'MinDomain')
    OPTIONS.MinDomain = -5 * ones(1, OPTIONS.Dim);
end
if ~isfield(OPTIONS, 'MaxDomain')
    OPTIONS.MaxDomain = +5 * ones(1, OPTIONS.Dim);
end
WorstRes = 0.25;
OPTIONS.BitsPerDim = ceil(log2((OPTIONS.MaxDomain - OPTIONS.MinDomain) / WorstRes + 1)); % bits per gene
OPTIONS.Resolution = (OPTIONS.MaxDomain - OPTIONS.MinDomain) ./ (2 .^ OPTIONS.BitsPerDim - 1);
OPTIONS.numVar = sum(OPTIONS.BitsPerDim); % number of bits in each chromosome
OPTIONS.OrderDependent = true;
% Initialize population
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    Population(popindex).chrom = round(rand(1, OPTIONS.numVar));
    if OPTIONS.Gray
        BitCountAccum = 0;
        for i = 1 : OPTIONS.Dim
            LowerBitIndex = BitCountAccum + 1;
            BitCountAccum = BitCountAccum + OPTIONS.BitsPerDim(i);
            UpperBitIndex = BitCountAccum;
            Bits = Population(popindex).chrom(LowerBitIndex : UpperBitIndex);
            Population(popindex).chrom(LowerBitIndex : UpperBitIndex) = BinaryToGray(Bits);
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
for popindex = 1 : length(Population)
    sum1 = 0;
    sum2 = 0;
    BitCountAccum = 0;
    for i = 1 : OPTIONS.Dim
        LowerBitIndex = BitCountAccum + 1;
        BitCountAccum = BitCountAccum + OPTIONS.BitsPerDim(i);
        UpperBitIndex = BitCountAccum;
        Bits = Population(popindex).chrom(LowerBitIndex : UpperBitIndex);
        gene = BitsToGene(Bits, OPTIONS.MinDomain(i), OPTIONS.Resolution(i), OPTIONS.Gray);
        sum1 = sum1 + gene^2;
        sum2 = sum2 + cos(2*pi*gene);
    end
    Population(popindex).cost = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/OPTIONS.Dim)) - exp(sum2/OPTIONS.Dim);    
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Feasible(Population, OPTIONS)
% Replace any genes outside of the allowable domain with randomly generated genes
for popindex = 1 : length(Population)
    BitCountAccum = 0;
    for i = 1 : OPTIONS.Dim
        LowerBitIndex = BitCountAccum + 1;
        BitCountAccum = BitCountAccum + OPTIONS.BitsPerDim(i);
        UpperBitIndex = BitCountAccum;
        while true
            Bits = Population(popindex).chrom(LowerBitIndex : UpperBitIndex);
            gene = BitsToGene(Bits, OPTIONS.MinDomain(i), OPTIONS.Resolution(i), OPTIONS.Gray);
            if (gene > OPTIONS.MinDomain(i)) && (gene < OPTIONS.MaxDomain(i)), break, end
            Population(popindex).chrom(LowerBitIndex : UpperBitIndex) = ...
                round(rand(1, OPTIONS.BitsPerDim(i))); % Replace gene with randomly generated gene
        end
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gene] = BitsToGene(Bits, MinDomain, Resolution, GrayFlag)
% Convert the bits in the Bits array to a real-valued gene
if GrayFlag
    Bits = GrayToBinary(Bits);
end
Gene = sum(Bits .* (2.^(0 : length(Bits)-1)));
Gene = MinDomain + (Gene - 1) * Resolution;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Binary] = GrayToBinary(Gray)
% Convert the gray-coded bit string in Gray to a binary-coded bit string in Binary.
% The first bit in each array is the LSB.
Binary = zeros(size(Gray));
Binary(end) = Gray(end);
for i = length(Gray)-1 : -1 : 1
    Binary(i) = mod(Binary(i+1) + Gray(i), 2);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gray] = BinaryToGray(Binary)
% Convert the binary-coded bit string in Binary to a gray-coded bit string in Gray.
% The first bit in each array is the LSB.
Gray = Binary;
for i = length(Gray)-1 : -1 : 1
    Gray(i) = xor(Binary(i+1), Binary(i));
end
return