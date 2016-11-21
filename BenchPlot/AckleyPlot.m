function AckleyPlot

xlimits = [-30 30];
dx = (xlimits(2) - xlimits(1)) / 250;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        sum1 = x^2 + y^2;
        sum2 = cos(2*pi*x) + cos(2*pi*y);
        obj(i,k) = 20 + exp(1) - 20 * exp(-0.2*sqrt(sum1/2)) - exp(sum2/2);
    end
end
close all;
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
xlabel('x:[-3 3]')
ylabel('y:[-3 3]')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])