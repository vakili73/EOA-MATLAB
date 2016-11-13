function [Population, NumDups] = ClearDups(Population, OPTIONS)
% Make sure there are no duplicate individuals in the population. This logic does not make 100% 
% sure that no duplicates exist, but any duplicates that are found are replaced with random individuals, 
% so there should be a good chance that there are no duplicates after this routine finishes.
NumDups = 0;
for i = 1 : length(Population)
    Chrom1 = Population(i).chrom;
    if ~OPTIONS.OrderDependent
        Chrom1 = sort(Chrom1);
    end
    for j = i+1 : length(Population)
        Chrom2 = Population(j).chrom;
        if ~OPTIONS.OrderDependent
            Chrom2 = sort(Chrom2);
        end
        if isequal(Chrom1, Chrom2)
            if isfield(OPTIONS, 'BitsPerDim')
                chrom = randi([0,1], [1,OPTIONS.numVar]);
            else
                chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
            end
            Population(j).chrom = chrom;
            Population(j) = OPTIONS.CostFunction(Population(j), OPTIONS);
            NumDups = NumDups + 1;
        end
    end
end
return