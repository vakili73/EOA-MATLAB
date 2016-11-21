function NPSOMonte

close all
nMonte = 20;
GenLimit = 50;
c4Array = 0.0 : 0.5 : 1.5+eps;
c5Array = 0.0 : 0.5 : 1.5+eps;
c6Array = 0.0 : 0.5 : 1.5+eps;
MinCostc4 = zeros(GenLimit+1, nMonte, length(c4Array));
MinCostc5 = zeros(GenLimit+1, nMonte, length(c5Array));
MinCostc6 = zeros(GenLimit+1, nMonte, length(c6Array));
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for j = 1 : length(c4Array)
        MinCostc4(:,i,j) = NPSO(@Schwefel226, false, 2, 2, 2, c4Array(j), 0, 0, 0.9, GenLimit); 
    end
    for j = 1 : length(c5Array)
        MinCostc5(:,i,j) = NPSO(@Schwefel226, false, 2, 2, 2, 0, c5Array(j), 0, 0.9, GenLimit); 
    end
    for j = 1 : length(c6Array)
        MinCostc6(:,i,j) = NPSO(@Schwefel226, false, 2, 2, 2, 0, 0, c6Array(j), 0.9, GenLimit); 
    end
end
colors = ['k: '; 'r- '; 'b--'; 'm-.'; 'ko '];
ncolors = size(colors, 1);
% Plot PSO performance for different c4 values
MinCostc4 = mean(MinCostc4, 2);
figure, hold on
Lgnd = cell(length(c4Array), 1);
for i = 1 : length(c4Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc4(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\phi_{4,\max} = ', num2str(c4Array(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
xlabel('Generation')
ylabel('Minimum Cost')
% Plot PSO performance for different c5 values
MinCostc5 = mean(MinCostc5, 2);
figure, hold on
Lgnd = cell(length(c5Array), 1);
for i = 1 : length(c5Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc5(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\phi_{5,\max} = ', num2str(c5Array(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
xlabel('Generation')
ylabel('Minimum Cost')
% Plot PSO performance for different c6 values
MinCostc6 = mean(MinCostc6, 2);
figure, hold on
Lgnd = cell(length(c6Array), 1);
for i = 1 : length(c6Array)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCostc6(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\phi_{6,\max} = ', num2str(c6Array(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
xlabel('Generation')
ylabel('Minimum Cost')