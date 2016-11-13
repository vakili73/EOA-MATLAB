function [ randVectorC ] = SingerPlot
% where a is a control parameter in (0.9, 1.08)
OPTIONS.noRandC = 10000;
x0 = @rand;
OPTIONS.x0 = feval(x0);
OPTIONS.a = linspace(0.9, 1.08, 100);

ia = 1;
while ia <= length(OPTIONS.a)
    a = OPTIONS.a(ia);
    randVectorC = zeros(1, OPTIONS.noRandC);
    for i = 1 : OPTIONS.noRandC
        randVectorC(i) = a * (7.86 * OPTIONS.x0 - 23.31 * (OPTIONS.x0 ^ 2) + ...
            28.75 * (OPTIONS.x0 ^ 3) - 13.302875 * (OPTIONS.x0 ^ 4));
        OPTIONS.x0 = feval(x0);
    end
    clf
    disp(['Constant a: ', num2str(a)]);
    subplot(1, 2, 1);
    scatter(1 : length(randVectorC), randVectorC, '.');
    title('Singer Chaotic Map Function');
    xlabel('Iteration');
    ylabel('Randomly Generated');
    
    subplot(1, 2, 2);
    histogram(randVectorC, 100);
    title('Singer Chaotic Map Function');
    xlabel('100 Bins');
    ylabel('Number Of Example');
    pause(0.1);
    ia = ia + 1;
end

end

