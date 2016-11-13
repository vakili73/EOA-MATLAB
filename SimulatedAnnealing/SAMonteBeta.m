function SAMonteBeta
nMonte = 20;
Problem = @Ackley;
Display = false;
GenLimit = 10000;
Restart = inf;
alpha = 0;
betaArr = [0.0002, 0.0010, 0.0005];
MinCostBeta1 = zeros(GenLimit+1, nMonte);
MinCostBeta2 = zeros(GenLimit+1, nMonte);
MinCostBeta3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostBeta1(:,i) = SA(Problem, Display, GenLimit, Restart, alpha, betaArr(1));
    MinCostBeta2(:,i) = SA(Problem, Display, GenLimit, Restart, alpha, betaArr(2));
    MinCostBeta3(:,i) = SA(Problem, Display, GenLimit, Restart, alpha, betaArr(3));
end
MinCostBeta1 = mean(MinCostBeta1, 2);
MinCostBeta2 = mean(MinCostBeta2, 2);
MinCostBeta3 = mean(MinCostBeta3, 2);
SetPlotOptions
figure
plot(0:GenLimit,MinCostBeta1,'r--', 0:GenLimit,MinCostBeta2,'k:', 0:GenLimit,MinCostBeta3,'b-');
xlabel('Iteration')
ylabel('Best Cost So Far')
legend(['beta = ', num2str(betaArr(1))], ['beta = ', num2str(betaArr(2))], ['beta = ', num2str(betaArr(3))]);