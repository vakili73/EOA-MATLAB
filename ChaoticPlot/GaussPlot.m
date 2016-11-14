function GaussPlot
% Gauss map function is defined by following equation
OPTIONS_.noRand = 100;

OPTIONS = struct();
OPTIONS.x0 = @rand;

randVector = Gauss(1, OPTIONS_.noRand, OPTIONS);

subplot(1, 2, 1);
scatter(1 : length(randVector), randVector, '*');
title('Gauss Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(randVector, 100);
title('Gauss Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

