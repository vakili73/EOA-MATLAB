function PSOFullyMonte
close all
nMonte = 20;
phiArray = [1, 2, 10, 20];
KArray = [0.1 0.5 0.9];
GenLimitPhi = 40;
MinCostphi = zeros(GenLimitPhi+1, nMonte, length(phiArray));
GenLimitK = 40;
MinCostK = zeros(GenLimitK+1, nMonte, length(KArray));
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for j = 1 : length(phiArray)
        MinCostphi(:,i,j) = PSOFully(@Ackley, false, phiArray(j), 0.9, GenLimitPhi);
    end
    for j = 1 : length(KArray)
        MinCostK(:,i,j) = PSOFully(@Ackley, false, 2, KArray(j), GenLimitK); 
    end
end
colors = ['k: '; 'r- '; 'b--'; 'm-.'; 'ko '];
ncolors = size(colors, 1);
% Plot PSO performance for different phiMax values
MinCostphi = mean(MinCostphi, 2);
figure, hold on
Lgnd = cell(length(phiArray), 1);
for i = 1 : length(phiArray)
    colorndx = mod(i, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimitPhi, MinCostphi(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\phi_{j,\max} = ', num2str(phiArray(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
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
    plot(0:GenLimitK, MinCostK(:, i), colors(colorndx, :))
    Lgnd{i} = ['$$\alpha = ', num2str(KArray(i)), '$$'];
end
h = legend(Lgnd);
set(h, 'Interpreter', 'latex')
xlabel('Generation')
ylabel('Minimum Cost')
box on