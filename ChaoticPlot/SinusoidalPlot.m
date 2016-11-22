function SinusoidalPlot
% where a is a control parameter
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(-5, 5, 1000);

OPTIONS = struct();
% OPTIONS.x0 = @rand;

for a = OPTIONS_.a
%     OPTIONS.a = a;
    randVector = Sinusoidal(1, OPTIONS_.noRand, OPTIONS)-Sinusoidal(1, OPTIONS_.noRand, OPTIONS);
    %% ploting
    clf
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('Sinusoidal Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVector, 100);
    title('Sinusoidal Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.01);
end

end

