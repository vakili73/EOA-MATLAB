function ICMICPlot
% where a is a control parameter
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(-1, 1, 100);

OPTIONS = struct();
OPTIONS.x0 = @randn;

for a = OPTIONS_.a
    OPTIONS.a = a;
    randVector = ICMIC(1, OPTIONS_.noRand, OPTIONS);
    %% ploting
    clf;
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('ICMIC Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVector, 100);
    title('ICMIC Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.1);
end

end

