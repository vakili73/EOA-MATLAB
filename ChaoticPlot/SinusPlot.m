function SinusPlot
% This chaotic function is defined by equation below
OPTIONS_.noRand = 10000;

OPTIONS = struct();
OPTIONS.x0 = @randn;

randVector = Sinus(1, OPTIONS_.noRand, OPTIONS);

subplot(1, 2, 1);
scatter(1 : length(randVector), randVector, '.');
title('Sinus Chaotic Map Function');
xlabel('Iteration');
ylabel('Randomly Generated');

subplot(1, 2, 2);
histogram(real(randVector), 100);
title('Sinus Chaotic Map Function');
xlabel('100 Bins');
ylabel('Number Of Example');

end

