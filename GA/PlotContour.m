function PlotContour(OPTIONS, Population, PlotNumber)
% Plots the population on top of the Ackley contour plot
if (PlotNumber == 1)
    figure
end
subplot(2, 2, PlotNumber);
[xArr, yArr, objArr] = AckleyContour(OPTIONS);
contour(xArr, yArr, objArr), hold on
gene = zeros(1, OPTIONS.Dim);
for popindex = 1 : OPTIONS.popsize
    BitCountAccum = 0;
    for i = 1 : OPTIONS.Dim
        LowerBitIndex = BitCountAccum + 1;
        BitCountAccum = BitCountAccum + OPTIONS.BitsPerDim(i);
        UpperBitIndex = BitCountAccum;
        Bits = Population(popindex).chrom(LowerBitIndex : UpperBitIndex);
        if OPTIONS.Gray
            Bits = GrayToBinary(Bits);
        end
        gene(i) = sum(Bits .* (2.^(0 : length(Bits)-1)));
        gene(i) = OPTIONS.MinDomain(i) + (gene(i) - 1) * OPTIONS.Resolution(i);
    end
    plot(gene(1), gene(2), '*')
end
if (PlotNumber == 1)
    xlabel('1st generation')
elseif (PlotNumber == 2)
    xlabel('4th generation')
elseif (PlotNumber == 3)
    xlabel('7th generation')
elseif (PlotNumber == 4)
    xlabel('10th generation')
end

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
