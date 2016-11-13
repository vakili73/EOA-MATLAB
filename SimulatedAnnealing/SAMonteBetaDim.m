function SAMonteBetaDim
nMonte = 20;
Problem = @AckleyScaled;
Display = false;
GenLimit = 10000;
betaArr = [0.005, 0.005;
    0.001, 0.005; 
    0.001, 0.001;
    0.005, 0.001];
Color = ['r--'; 'b- '; 'k: '; 'm-.'];
MinCostBeta = zeros(size(betaArr, 1), GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for k = 1 : size(betaArr, 1)
        MinCostBeta(k, :, i) = SADimension(Problem, Display, GenLimit, betaArr(k, 1), betaArr(k, 2));
    end
end
MinCostBeta = mean(MinCostBeta, 3);
SetPlotOptions
figure, hold on
for k = 1 : size(betaArr, 1)
    plot(0:GenLimit, MinCostBeta(k,:), Color(k,:))
end
xlabel('Iteration')
ylabel('Best Cost So Far')
legend('beta = 0.005', 'beta = 0.001 (odd), 0.005 (even)', 'beta = 0.001', 'beta = 0.005 (odd), 0.001 (even)');
