function LogisticPlot
% where a is a control parameter in (0, 1)
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(-10, 10, 100);

OPTIONS = struct();
% OPTIONS.x0 = @rand;

for a = OPTIONS_.a
%     OPTIONS.a = a;
    randVector = Logistic(1, OPTIONS_.noRand, OPTIONS)-Logistic(1, OPTIONS_.noRand, OPTIONS);
    %% ploting
    clf
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('Logistic Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVector, 100);
    title('Logistic Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.1);
end

end

