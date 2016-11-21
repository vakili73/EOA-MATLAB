function ChebyshevPlot
% where a is a control parameter in (0, 1)
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(-10, 10, 100);

OPTIONS = struct();
% OPTIONS.x0 = @rand;

for a = OPTIONS_.a
    OPTIONS.a = a;
    randVector = Chebyshev(1, OPTIONS_.noRand, OPTIONS)-Chebyshev(1, OPTIONS_.noRand, OPTIONS);
%     randVector = mapminmax(randVector);
    %% ploting
    clf
    disp(['Constant a: ', num2str(a), ' Min: ', num2str(min(randVector)), ' Max: ', num2str(max(randVector)) ]);    
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('Chebyshev Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(real(randVector), 100);
    title('Chebyshev Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.5);
end

end

