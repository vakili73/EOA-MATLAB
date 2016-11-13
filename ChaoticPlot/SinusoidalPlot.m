function [ randVectorC ] = SinusoidalPlot
% where a is a control parameter
OPTIONS.noRandC = 10000;
x0 = @rand;
OPTIONS.x0 = feval(x0);
OPTIONS.a = linspace(-1, 1, 100);

ia = 1;
while ia <= length(OPTIONS.a)
    a = OPTIONS.a(ia);
    randVectorC = zeros(1, OPTIONS.noRandC);
    for i = 1 : OPTIONS.noRandC
        randVectorC(i) = a * (OPTIONS.x0 ^ 2) * sin(pi * OPTIONS.x0);
        OPTIONS.x0 = feval(x0);
    end
    clf
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVectorC), randVectorC, '.');
    title('Sinusoidal Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVectorC, 100);
    title('Sinusoidal Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.1);
    ia = ia + 1;
end

end

