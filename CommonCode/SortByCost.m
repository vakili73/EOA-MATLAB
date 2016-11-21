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
