function [ randVectorC ] = TentPlot
OPTIONS.noRandC = 10000;
x0 = @randn;
OPTIONS.x0 = feval(x0);

randVectorC = zeros(1, OPTIONS.noRandC);
for i = 1 : OPTIONS.noRandC
    if OPTIONS.x0 < 0.5
        randVectorC(i) = 2 * OPTIONS.x0;
    end
    if OPTIONS.x0 >= 0.5
        randVectorC(i) = 2 * (1 - OPTIONS.x0);
    end
    OPTIONS.x0 = feval(x0);
end
subplot(1, 2, 1);
scatter(1 : length(randVectorC), randVectorC, '.');
title('Tent Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(randVectorC, 100);
title('Tent Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

