function [ randVectorC ] = SinusPlot
% This chaotic function is defined by equation below
OPTIONS.noRandC = 10000;
x0 = @rand;
OPTIONS.x0 = feval(x0);

randVectorC = zeros(1, OPTIONS.noRandC);
for i = 1 : OPTIONS.noRandC
    randVectorC(i) = 2.3 * (OPTIONS.x0 ^ (2 * sin(pi * OPTIONS.x0)));
    OPTIONS.x0 = feval(x0);
end
subplot(1, 2, 1);
scatter(1 : length(randVectorC), randVectorC, '.');
title('Sinus Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(randVectorC, 100);
title('Sinus Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

