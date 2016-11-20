function RandNormPlot
% This chaotic function is defined by equation below
OPTIONS_.noRand = 10000;
randVector = randn(1, OPTIONS_.noRand);
randVector = mapminmax(randVector);

subplot(1, 2, 1);
scatter(1 : length(randVector), randVector, '.');
title('RandNorm Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(real(randVector), 100);
title('RandNorm Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

disp(['Mean: ', num2str(mean(randVector)), '   STD: ', num2str(std(randVector))]);
end

