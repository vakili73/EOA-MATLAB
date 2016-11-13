function ShekelPlot

xlimits = [-50 50];
dx = (xlimits(2) - xlimits(1)) / 100;
N = floor((xlimits(2) - xlimits(1)) / dx) + 1;
a = zeros(2, 25);
b0 = [-32 -16 0 16 32];
a(1, :) = [b0 b0 b0 b0 b0];
b = -32 * [1 1 1 1 1];
a(2, :) = [b b+16 b+32 b+48 b+64];
obj = zeros(N, N);
i = 0;
for x = xlimits(1) : dx : xlimits(2)
    i = i + 1;
    k = 0;
    for y = xlimits(1) : dx : xlimits(2)
        k = k + 1;
        for p = 1 : size(a, 2)
            obj(i,k) = obj(i,k) + 1 / (p + (x - a(1,p))^6 + (y - a(2,p))^6);
        end
        obj(i,k) = 1 / (1/500 + obj(i,k));
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