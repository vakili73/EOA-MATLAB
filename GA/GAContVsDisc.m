function GAContVsDisc
% Compare the performance of the continuous GA to that of the discrete GA
nMonte = 50;
Display = false;
GenLimit = 20;
Gray = false;
NumElites = 1;
PopSize = 10;
PmutateDisc = 0.02;
PmutateCont = 4 * PmutateDisc; % 4 bits per dimension in discrete problem (depends on search domain)
MinDomain = -1;
MaxDomain = +1;
Dimensions = 2;
ContCostArr = zeros(GenLimit+1, nMonte);
DiscCostArr = zeros(GenLimit+1, nMonte);
ContAvgCostArr = zeros(GenLimit+1, nMonte);
DiscAvgCostArr = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    RandSeed = fix(sum(100*clock));
    [DiscCostArr(:,i), DiscAvgCostArr(:,i)] = GA(@AckleyDisc, Display, RandSeed, GenLimit, Gray, NumElites, ...
        PopSize, PmutateDisc, MinDomain, MaxDomain, Dimensions);
    [ContCostArr(:,i), ContAvgCostArr(:,i)] = GA(@Ackley, Display, RandSeed, GenLimit, Gray, NumElites, ...
        PopSize, PmutateCont, MinDomain, MaxDomain, Dimensions);
end
DiscCostArr = mean(DiscCostArr, 2);
ContCostArr = mean(ContCostArr, 2);
DiscAvgCostArr = mean(DiscAvgCostArr, 2);
ContAvgCostArr = mean(ContAvgCostArr, 2);
% Plot minimum cost
figure, hold on, box on
plot(0:GenLimit, DiscCostArr, 'r--')
plot(0:GenLimit, ContCostArr, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Binary GA', 'Continuous GA')
% Plot average cost
figure, hold on, box on
plot(0:GenLimit, DiscAvgCostArr, 'r--')
plot(0:GenLimit, ContAvgCostArr, 'b-')
xlabel('Generation')
ylabel('Average Cost')
legend('Binary GA', 'Continuous GA')