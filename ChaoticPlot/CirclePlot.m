function [ randVectorC ] = CirclePlot
% By x0 = 0.5 and b = 0.2, chaotic numbers will be generated in (0,1)
OPTIONS.noRandC = 10000;
x0 = @randn;
OPTIONS.x0 = feval(x0);
OPTIONS.a = linspace(-100, 100, 1000);
OPTIONS.b = linspace(0, 1, 1);

ia = 1;
while ia <= length(OPTIONS.a)
    a = OPTIONS.a(ia);
    ib = 1;
    while ib <= length(OPTIONS.b)
        b = OPTIONS.b(ib);
        randVectorC = zeros(1, OPTIONS.noRandC);
        for i = 1 : OPTIONS.noRandC
            randVectorC(i) = OPTIONS.x0 + b - a / (2 * pi) * sin(2 * pi * OPTIONS.x0);
            OPTIONS.x0 = feval(x0);
        end
        clf
        disp(['Constant a: ', num2str(a), ' --- Constant b:', num2str(b)]);
        subplot(1, 2, 1);
        scatter(1 : length(randVectorC), randVectorC, '.');
        title('Circle Chaotic Map Function');
        xlabel('Iteration');
        ylabel('Randomly Generated');
        
        subplot(1, 2, 2);
        histogram(randVectorC, 100);
        title('Circle Chaotic Map Function');
        xlabel('100 Bins');
        ylabel('Number Of Example');
        pause(0.01);
        ib = ib + 1;
    end
    ia = ia + 1;
end

end

