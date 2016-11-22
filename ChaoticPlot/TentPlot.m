function TentPlot
OPTIONS_.noRand = 10000;

OPTIONS = struct();
% OPTIONS.x0 = @rand;

randVector = Tent(1, OPTIONS_.noRand, OPTIONS)-Tent(1, OPTIONS_.noRand, OPTIONS);

subplot(1, 2, 1);
scatter(1 : length(randVector), randVector, '.');
title('Tent Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(randVector, 100);
title('Tent Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

