function CauchyGaussian

% Plot a Cauchy PDF and a Gaussian PDF

x = -5 : .01 : 5;
sigma = 1;
gamma = sigma * sqrt(2 / pi);

Cauchy = gamma ./ (x.^2 + gamma^2) / pi;
Gauss = exp(-x.^2 / 2 / sigma.^2) / sigma / sqrt(2 * pi);

close all
figure
plot(x,Gauss,'r-', x,Cauchy,'b--')
xlabel('x')
ylabel('PDF(x)')
legend('Gaussian PDF', 'Cauchy PDF')