function QuarticPlot

xlimits = [-1.28 1.28];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
obj = zeros(N, N);
i = 0;
r = rand;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        obj(i,k) = x^4 + 2*y^4 + r;
    end
end
close all;
[x, y] = meshgrid(xlimits(1) : dx : xlimits(2));
meshc(x,y,obj);
zlabel('f(x)')
set(gca, 'XTick', [])
set(gca, 'YTick', [])
set(gca, 'ZTick', [])
temp = axis;
axis([xlimits, xlimits, temp(5), max(max(obj))])