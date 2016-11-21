function PSOMonte

close all
nMonte = 20;
GenLimit = 50;
c1Array = 1 : .5 : 3+eps;
c2Array = 1 : .5 : 3+eps;
c3Array = 1 : .5 : 3+eps;
KArray = 0.1 : 0.4 : 1.3+eps;
MinCostc1 = zeros(GenLimit+1, nMonte, length(c1Array));
MinCostc2 = zeros(GenLimit+1, nMonte, length(c2Array));
MinCostc3 = zeros(GenLimit+1, nMonte, length(c3Array));
MinCostK = zeros(GenLimit+1, nMonte, length(KArray));
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for j = 1 : length(c1Array)
        MinCostc1(:,i,j) = PSO(@Ackley, false, c1Array(j), 2.1, 2.1, 0.9, GenLimit); 
    end
    for j = 1 : length(c2Array)
        MinCostc2(:,i,j) = PSO(@Ackley, false, 2.1, c2Array(j), 2.1, 0.9, GenLimit); 
    end
    for j = 1 : length(c3Array)
        MinCostc3(:,i,j) = PSO(@Ackley, false, 2.1, 2.1, c3Array(j), 0.9, GenLimit); 
    end
    for j = 1 : length(KArray)
        MinCostK(:,i,j) = PSO(@Ackley, false, 2.1, 2.1, 2.1, KArray(j), GenLimit); 
    end
end
colors = ['k: '; 'r- '; 'b--'; 'm-.'; 'ko '];
ncolors = size(colors, 1);
% Plot PSO performance for different c1 values
MinCostc1 = mean(MinCostc1, 2);
figure, hold on
Lgnd = cell(length(c1Array), 1);
for i = 1 : length(c1Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc1(:, i), colors(colorndx, :))
    Lgnd{i} = ['c_1 = ', num2str(c1Array(i))];
end
legend(Lgnd);
xlabel('Generation')
ylabel('Minimum Cost')
% Plot PSO performance for different c2 values
MinCostc2 = mean(MinCostc2, 2);
figure, hold on
Lgnd = cell(length(c2Array), 1);
for i = 1 : length(c2Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc2(:, i), colors(colorndx, :))
    Lgnd{i} = ['c_2 = ', num2str(c2Array(i))];
end
legend(Lgnd);
xlabel('Generation')
ylabel('Minimum Cost')
% Plot PSO performance for different c3 values
MinCostc3 = mean(MinCostc3, 2);
figure, hold on
Lgnd = cell(length(c3Array), 1);
for i = 1 : length(c3Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc3(:, i), colors(colorndx, :))
    Lgnd{i} = ['c_3 = ', num2str(c3Array(i))];
end
legend(Lgnd);
xlabel('Generation')
ylabel('Minimum Cost')
% Plot PSO performance for different K values
MinCostK = mean(MinCostK, 2);
figure, hold on
Lgnd = cell(length(KArray), 1);
for i = 1 : length(KArray)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostK(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\alpha = ', num2str(KArray(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
xlabel('Generation')
ylabel('Minimum Cost')