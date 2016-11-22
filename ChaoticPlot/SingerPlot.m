function SingerPlot
% where a is a control parameter in (0.9, 1.08)
OPTIONS_.noRand = 10000;
OPTIONS_.a = linspace(0.9, 1.08, 100);

OPTIONS = struct();
% OPTIONS.x0 = @rand;

for a = OPTIONS_.a
%     OPTIONS.a = a;
    randVector = Singer(1, OPTIONS_.noRand, OPTIONS)-Singer(1, OPTIONS_.noRand, OPTIONS);
    %% ploting
    clf
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVector), randVector, '.');
    title('Singer Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVector, 100);
    title('Singer Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.1);
end

end

