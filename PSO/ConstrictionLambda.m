function ConstrictionLambda(phi)

if ~exist('phi', 'var')
    phi = 5;
end

alphaBound = 2 / (phi - 2);
alpha = 0 : alphaBound/1000 : alphaBound;
lambda1 = (1/2) * (1 - alpha * (phi - 1) + sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
lambda2 = (1/2) * (1 - alpha * (phi - 1) - sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
close all
figure, hold on
plot(alpha, abs(lambda1), 'r-', alpha, abs(lambda2), 'b--')
xlabel('$$K$$', 'Interpreter', 'latex')
ylabel('Magnitude')
h = legend('$$|\lambda_1|$$', '$$|\lambda_2|$$');
set(h, 'Interpreter', 'latex')
% Finish the lambda1 part of the plot
temp = axis;
alpha = alphaBound : (temp(2)-alphaBound)/1000 : temp(2);
lambda1 = (1/2) * (1 - alpha * (phi - 1) + sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
plot(alpha, abs(lambda1), 'r-')

K1 = (phi + 1 - 2 * sqrt(phi)) / (phi  - 1)^2;
lambda1 = (sqrt(phi) - 1) / (phi - 1);
arrow([.2, .2], [K1, lambda1], 'Width', 2) 
str(1) = {'$$K=K_1$$'};
str(2) = {'$$|\lambda| = (\sqrt{\phi}-1)/(\phi-1)$$'};
text('Position', [.2, .2], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', str);

K2 = (phi + 1 + 2 * sqrt(phi)) / (phi  - 1)^2;
lambda1 = (1 + sqrt(phi)) / (phi - 1);
arrow([.45, .9], [K2, lambda1], 'Width', 2)
str(1) = {'$$K=K_2$$'};
str(2) = {'$$|\lambda| = (\sqrt{\phi}+1)/(\phi-1)$$'};
text('Position', [.1, .9], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', str);

K3 = 2 / (phi - 2);
lambda2 = 1;
arrow([.7, .8], [K3, lambda2], 'Width', 2)
str(1) = {'$$K=K_3$$'};
str(2) = {'$$\lambda_2 = -1$$'};
text('Position', [.7, .8], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', str);

% Next plot

alphaBound = 2 * (phi + 1 + 2 * sqrt(phi)) / ((phi - 1)^2);
alpha = 0 : alphaBound/100000 : alphaBound;
lambda1 = (1/2) * (1 - alpha * (phi - 1) + sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
lambda2 = (1/2) * (1 - alpha * (phi - 1) - sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
figure; hold on; grid
plot(real(lambda1), imag(lambda1), 'r-')
plot(real(lambda2), imag(lambda2), 'b--')
h = legend('$$\lambda_1$$', '$$\lambda_2$$');
set(h, 'Interpreter', 'latex')
xlabel('$$\mbox{Real}(\lambda)$$', 'Interpreter', 'latex')
ylabel('$$\mbox{Imag}(\lambda)$$', 'Interpreter', 'latex')
lambda1 = lambda1(end) : 0.01 : 1/(1-phi);
plot(lambda1, zeros(size(lambda1)), 'r-') % finish the lambda1 plot to its limit
% Plot the unit circle
x = -1 : .01 : 1;
y = sqrt(1 - x.^2);
plot(x, y, 'k:', x, -y, 'k:')
text(-1.2, 0.9, 'Unit Circle', 'BackgroundColor', 'w')
axis([-1.4 1.4 -1.2 1.2])
arrow([1, -.7], [1, 0], 'Width', 2) % arrow to lambda1(K=0) = 1
text('Position', [.9, -.7], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda_1(0)$$')
arrow([0, -.7], [0, 0], 'Width', 2) % arrow to lambda2(K=1) = 0
text('Position', [-.1, -.7], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda_2(0)$$')
lambda1 = (1 - sqrt(phi)) / (1 - phi);
arrow([lambda1, -.7], [lambda1, 0], 'Width', 2) % arrow to lambda1=lambda2>0
text('Position', [lambda1-.1, -.7], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda(K_1)$$')
lambda1 = (1 + sqrt(phi)) / (1 - phi);
arrow([lambda1, -.7], [lambda1, 0], 'Width', 2) % arrow to lambda1=lambda2<0
text('Position', [lambda1-.05, -.7], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda(K_2)$$')
arrow([1/(1-phi), -.3], [1/(1-phi), 0], 'Width', 2) % arrow to lambda1(K=inf)
text('Position', [-.5, -.35], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda_1(\infty)$$')
arrow([-1, -.7], [-1, 0], 'Width', 2) % arrow to lambda2=-1
text('Position', [-1.3, -.7], 'BackgroundColor', 'w', 'Interpreter', 'latex', 'String', '$$\lambda_2(K_3)$$')
grid off

return

syms phi alpha kappa real;
%alpha = 2 * kappa / (phi - 2 + sqrt(phi^2 - 4 * phi));
lambda1 = (1/2) * (1 - alpha * (phi - 1) + sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
lambda2 = (1/2) * (1 - alpha * (phi - 1) - sqrt(1 + alpha.^2 * (phi - 1)^2 - 2 * alpha * (phi + 1)));
pretty(lambda1)
pretty(lambda2)
pretty(diff(lambda1, 'alpha'))
pretty(diff(lambda2, 'alpha'))

return

lambda1Calc = matlabFunction(lambda1);
lambda2Calc = matlabFunction(lambda2);

clear phi kappa lambda1 lambda2
i = 0;
for phi = 4 : .01 : 6
    i = i + 1;
    k = 0;
    for kappa = 0 : .01 : 1
        k = k + 1;
        lambda1(i,k) = abs(lambda1Calc(kappa, phi));
        lambda2(i,k) = abs(lambda2Calc(kappa, phi));
    end
end
figure;
[phi, kappa] = meshgrid((0:.01:1), (4:.01:6));
mesh(kappa, phi, lambda1);
figure;
mesh(kappa, phi, lambda2);

phi = 4 : .01 : 6;
diff1 = 2 ./ (phi - 2 + sqrt(phi.^2 - 4 * phi)) - (phi + 1 - 2 * sqrt(phi)) ./ ((phi - 1).^2);
diff2 = (phi + 1 + 2 * sqrt(phi)) ./ ((phi - 1).^2) - 2 ./ (phi - 2 + sqrt(phi.^2 - 4 * phi));
figure
plot(phi, diff1, 'r', phi, diff2, 'b')
