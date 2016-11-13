function [ randVectorC ] = GaussPlot
% Gauss map function is defined by following equation
OPTIONS.noRandC = 100;
x0 = @rand;
OPTIONS.x0 = feval(x0);

randVectorC = zeros(1, OPTIONS.noRandC);
for i = 1 : OPTIONS.noRandC
    if OPTIONS.x0 == 0
        randVectorC(i) = 0;
    else
        randVectorC(i) = 1 / OPTIONS.x0;
    end
    OPTIONS.x0 = feval(x0);
end        
subplot(1, 2, 1);
scatter(1 : length(randVectorC), randVectorC, '*');
title('Gauss Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(randVectorC, 100);
title('Gauss Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

