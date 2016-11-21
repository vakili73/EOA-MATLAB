function Population = PopPick(Population, popsize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Temp = struct('chrom',[],'cost',[]);
for i = 1:popsize
    Temp(i) = Population(i);
end
Population = Temp;

end

