function BernoulliPlot
% where a is a control parameter in (0, 0.5)
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(0.01, 0.5, 100);

OPTIONS = struct();
% OPTIONS.x0 = @rand;

for a = OPTIONS_.a
%     OPTIONS.a = a;
%     randVector = Bernoulli(1, OPTIONS_.noRand, OPTIONS);
    randVector = Bernoulli(1, OPTIONS_.noRand, OPTIONS)-Bernoulli(1, OPTIONS_.noRand, OPTIONS);
%     randVector = rand(1, OPTIONS_.noRand)-rand(1, OPTIONS_.noRand);
%     randVector = randn(1, OPTIONS_.noRand);
    
    %% ploting
    clf;
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('Bernoulli Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVector, 100);
    title('Bernoulli Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.01);
end

end

