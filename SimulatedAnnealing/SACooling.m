function SACooling

SetPlotOptions

% Inverse cooling schedule
beta = [0.001, 0.005, 0.010];
kArr = 1 : 5000;
TArr = ones(length(beta), length(kArr));
for k = 2 : length(kArr)
    TArr(:, k) = TArr(:, k-1) ./ (1 + beta(:) .* TArr(:, k-1));
end
figure
plot(kArr,TArr(1,:),'r-', kArr,TArr(2,:),'b--', kArr,TArr(3,:),'k:')
xlabel('Iteration')
ylabel('Normalized Temperature')
legend(['beta = ', num2str(beta(1))], ['beta = ', num2str(beta(2))], ['beta = ', num2str(beta(3))])

% Exponential cooling schedule
a = [0.9996, 0.9994, 0.9992];
kArr = 1 : 5000;
TArr = ones(length(beta), length(kArr));
for k = 2 : length(kArr)
    TArr(:, k) = a(:) .* TArr(:, k-1);
end
figure
plot(kArr,TArr(1,:),'r-', kArr,TArr(2,:),'b--', kArr,TArr(3,:),'k:')
xlabel('Iteration')
ylabel('Normalized Temperature')
legend(['a = ', num2str(a(1))], ['a = ', num2str(a(2))], ['a = ', num2str(a(3))])

% Time-varying exponential cooling schedule
kArr = 1 : 50;
TArr = ones(1, length(kArr));
for k = 2 : length(kArr)
    TArr(:, k) = TArr(:, k-1) / k;
end
figure
plot(kArr,TArr(1,:),'r-')
xlabel('Iteration')
ylabel('Normalized Temperature')

% Logarithmic cooling schedule
kArr = 1 : 5000;
TArr = ones(1, length(kArr));
for k = 2 : length(kArr)
    TArr(:, k) = TArr(:, k) / log(k);
end
figure
subplot(2,1,1)
plot(kArr(1:50),TArr(1,1:50),'r-')
subplot(2,1,2)
plot(kArr,TArr(1,:),'r-')
xlabel('Iteration')
ylabel('Normalized Temperature')

